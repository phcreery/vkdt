/*
 * Image alignment module.
 * This module aligns a target image to a reference image.
 * It detects SURF features on the GPU.
 * The detector is a port of gpusurf by Paul Furgale and Chi Hay Tong.
 * It matches the features, and fits an affine matrix or a homography with RANSAC.
 * It applies the transformation to the target image in a compute shader.
 */

#include "modules/api.h"

#include <math.h>
#include <string.h>

#define SURF_OCTAVES   4
#define SURF_INTERVALS 4
#define CAND_MAX       16384
#define KP_STRIDE      8
#define DESC_DIM       64
#define FEAT_STRIDE    (KP_STRIDE + DESC_DIM)
#define MAX_FEAT       1024

// The number of threads for a one dimensional shader with 'ls' threads in each work group.
#define WG(n, ls) ((((n) + (ls) - 1) / (ls)) * DT_LOCAL_SIZE_X)

// Uniform parameters. The layout must match params_t in surf.glsl.
typedef struct shader_params_t
{
  float thresh;
  float ratio;
  float rthresh;
  int   maxfeat;
  int   model;
} shader_params_t;

// Push constants for the feature detection. The layout must match push_t in surf.glsl.
typedef struct geo_t
{
  int   size[SURF_OCTAVES][4];  // x size, y size, border, step
  float scale[SURF_OCTAVES][SURF_INTERVALS];
  int   offset[4];              // start of each octave in the hessian buffer
  int   dim[4];                 // integral image width, height, reduction factor
} geo_t;

// Push constants of the solver. The layout must match push_t in solve.comp.
typedef struct solve_pc_t
{
  int   dim[2];
  float factor[2];
} solve_pc_t;

// Compute the layout of the octaves for an image of the reduced size wd x ht.
// This follows GpuSurfOctave::init() of gpusurf with the default parameters.
// Return the number of values in the hessian buffer.
static int
compute_geo(
    geo_t *g,
    const int wd,
    const int ht,
    const int red)
{
  const float l1 = 3.0f/1.5f, l2 = 5.0f/1.5f, l3 = 3.0f/1.5f, l4 = 1.0f/1.5f;
  const float mask_width  = l2 * 0.5f;
  const float mask_height = 1.0f + l1;
  const float base_scale  = 2.0f;
  memset(g, 0, sizeof(*g));
  g->dim[0] = wd;
  g->dim[1] = ht;
  g->dim[2] = red;
  int total = 0;
  for(int o = 0; o < SURF_OCTAVES; o++)
  {
    const int step = 1 << o;
    const float d = base_scale * (1 << o) / (SURF_INTERVALS - 2);
    for(int i = 0; i < SURF_INTERVALS; i++)
      g->scale[o][i] = base_scale * (1 << o) + d * (i - 1.0f) + 0.5f;
    const float smax = g->scale[o][SURF_INTERVALS-1];
    const int border = (int)ceilf(smax * fmaxf(fmaxf(mask_width, mask_height), l3 + l4 * 0.5f));
    int xs = (wd - 2*border) / step;
    int ys = (ht - 2*border) / step;
    if(xs <= 0 || ys <= 0) xs = ys = 0; // the image is too small for this octave
    g->size[o][0] = xs;
    g->size[o][1] = ys;
    g->size[o][2] = border;
    g->size[o][3] = step;
    g->offset[o]  = total;
    total += xs * ys * SURF_INTERVALS;
  }
  return total;
}

// The factor by which the module reduces the image before the feature detection.
static int
reduction(
    const dt_roi_t *roi,
    const int       maxdim)
{
  const int m = roi->wd > roi->ht ? roi->wd : roi->ht;
  const int r = (m + maxdim - 1) / maxdim;
  return r < 1 ? 1 : r;
}

static int
param_int(
    dt_module_t *module,
    const char  *name)
{
  return dt_module_param_int(module, dt_module_get_param(module->so, dt_token(name)))[0];
}

static float
param_float(
    dt_module_t *module,
    const char  *name)
{
  return dt_module_param_float(module, dt_module_get_param(module->so, dt_token(name)))[0];
}

void
modify_roi_out(
    dt_graph_t  *graph,
    dt_module_t *module)
{
  module->connector[2].roi = module->connector[0].roi;
}

dt_graph_run_t
check_params(
    dt_module_t *module,
    uint32_t     parid,
    uint32_t     num,
    void        *oldval)
{
  const dt_token_t name = module->so->param[parid]->name;
  // these parameters change the size of the buffers
  if(name == dt_token("maxfeat") || name == dt_token("maxdim"))
    return s_graph_run_all;
  return s_graph_run_record_cmd_buf; // minimal parameter upload to uniforms
}

void
commit_params(
    dt_graph_t  *graph,
    dt_module_t *module)
{
  shader_params_t *p = (shader_params_t *)module->committed_param;
  p->thresh  = param_float(module, "thresh") * 1e-4f; // the user interface uses units of 1e-4
  p->ratio   = param_float(module, "ratio");
  p->rthresh = param_float(module, "rthresh");
  p->maxfeat = CLAMP(param_int(module, "maxfeat"), 8, MAX_FEAT);
  p->model   = param_int(module, "model");
}

int
init(dt_module_t *mod)
{
  mod->committed_param_size = sizeof(shader_params_t);
  return 0;
}

// Add the nodes that compute the SURF features of one input image.
// Return the id of the last node. Its connector 2 is the feature buffer.
// The feature buffer has one uint count, followed by the records of the features.
// Each record has the keypoint (x, y, size, response, angle, octave, laplacian sign)
// and the descriptor with 64 floats.
static int
add_surf_nodes(
    dt_graph_t  *graph,
    dt_module_t *module,
    const int    cid,      // input connector of the module
    const geo_t *geo,
    const int    nhess,    // number of values in the hessian buffer
    const int    maxfeat)
{
  const int wd = geo->dim[0];
  const int ht = geo->dim[1];
  const int pc = sizeof(geo_t);
  const int *pcd = (const int *)geo;
  const dt_roi_t ii_roi   = {.wd = (wd + 1) * (ht + 1), .ht = 1};
  const dt_roi_t hess_roi = {.wd = nhess, .ht = 1};
  const dt_roi_t cand_roi = {.wd = 1 + CAND_MAX * KP_STRIDE, .ht = 1};
  const dt_roi_t sel_roi  = {.wd = 1 + maxfeat * KP_STRIDE, .ht = 1};
  const dt_roi_t feat_roi = {.wd = 1 + maxfeat * FEAT_STRIDE, .ht = 1};

  const int id_row = dt_node_add(graph, module, "algnsurf", "iirow",
      WG(ht, 64), DT_LOCAL_SIZE_Y, 1, pc, pcd, 2,
      "input", "read",  "rgba", "*",   dt_no_roi,
      "rows",  "write", "ssbo", "u32", &ii_roi);
  const int id_col = dt_node_add(graph, module, "algnsurf", "iicol",
      WG(wd + 1, 64), DT_LOCAL_SIZE_Y, 1, pc, pcd, 2,
      "rows", "read",  "ssbo", "u32", dt_no_roi,
      "ii",   "write", "ssbo", "u32", &ii_roi);
  const int id_hess = dt_node_add(graph, module, "algnsurf", "hessian",
      geo->size[0][0], geo->size[0][1], SURF_OCTAVES * SURF_INTERVALS, pc, pcd, 2,
      "ii",   "read",  "ssbo", "u32", dt_no_roi,
      "hess", "write", "ssbo", "f32", &hess_roi);
  const int id_nms = dt_node_add(graph, module, "algnsurf", "nms",
      geo->size[0][0], geo->size[0][1], SURF_OCTAVES * (SURF_INTERVALS - 2), pc, pcd, 2,
      "hess", "read",  "ssbo", "f32", dt_no_roi,
      "cand", "write", "ssbo", "f32", &cand_roi);
  const int id_sel = dt_node_add(graph, module, "algnsurf", "select",
      128, 128, 1, 0, 0, 2,
      "cand", "read",  "ssbo", "f32", dt_no_roi,
      "sel",  "write", "ssbo", "f32", &sel_roi);
  const int id_desc = dt_node_add(graph, module, "algnsurf", "describe",
      maxfeat * DT_LOCAL_SIZE_X, DT_LOCAL_SIZE_Y, 1, pc, pcd, 3,
      "ii",   "read",  "ssbo", "u32", dt_no_roi,
      "sel",  "read",  "ssbo", "f32", dt_no_roi,
      "feat", "write", "ssbo", "f32", &feat_roi);

  dt_connector_copy(graph, module, cid, id_row, 0);
  CONN(dt_node_connect(graph, id_row,  1, id_col,  0));
  CONN(dt_node_connect(graph, id_col,  1, id_hess, 0));
  CONN(dt_node_connect(graph, id_col,  1, id_desc, 0));
  CONN(dt_node_connect(graph, id_hess, 1, id_nms,  0));
  CONN(dt_node_connect(graph, id_nms,  1, id_sel,  0));
  CONN(dt_node_connect(graph, id_sel,  1, id_desc, 1));
  graph->node[id_nms].connector[1].flags |= s_conn_clear; // the candidate counter starts at zero
  return id_desc;
}

void
create_nodes(
    dt_graph_t  *graph,
    dt_module_t *module)
{
  if(!dt_connected(module->connector + 1))
  { // there is nothing to align
    dt_connector_bypass(graph, module, 0, 2);
    return;
  }
  const dt_roi_t roi_ref = module->connector[0].roi;
  const dt_roi_t roi_tgt = module->connector[1].roi;
  const dt_roi_t roi_out = module->connector[2].roi;
  const int maxfeat = CLAMP(param_int(module, "maxfeat"), 8, MAX_FEAT);
  const int maxdim  = param_int(module, "maxdim");

  const int red_ref = reduction(&roi_ref, maxdim);
  const int red_tgt = reduction(&roi_tgt, maxdim);
  geo_t geo_ref, geo_tgt;
  const int nhess_ref = compute_geo(&geo_ref, roi_ref.wd / red_ref, roi_ref.ht / red_ref, red_ref);
  const int nhess_tgt = compute_geo(&geo_tgt, roi_tgt.wd / red_tgt, roi_tgt.ht / red_tgt, red_tgt);
  if(geo_ref.size[0][0] < 3 || geo_tgt.size[0][0] < 3)
  { // the image is too small for the largest filter
    snprintf(module->graph->gui_msg_buf, sizeof(module->graph->gui_msg_buf),
        "algnsurf: the images are too small");
    module->graph->gui_msg = module->graph->gui_msg_buf;
    dt_connector_bypass(graph, module, 1, 2);
    return;
  }

  const dt_roi_t mt_roi = {.wd = maxfeat, .ht = 1};
  const dt_roi_t h_roi  = {.wd = 13, .ht = 1};

  const int id_ref = add_surf_nodes(graph, module, 0, &geo_ref, nhess_ref, maxfeat);
  const int id_tgt = add_surf_nodes(graph, module, 1, &geo_tgt, nhess_tgt, maxfeat);

  const int id_match = dt_node_add(graph, module, "algnsurf", "match",
      WG(maxfeat, 64), DT_LOCAL_SIZE_Y, 1, 0, 0, 3,
      "ref",   "read",  "ssbo", "f32", dt_no_roi,
      "tgt",   "read",  "ssbo", "f32", dt_no_roi,
      "match", "write", "ssbo", "u32", &mt_roi);
  const solve_pc_t spc = {
    .dim    = {geo_ref.dim[0], geo_ref.dim[1]},
    .factor = {(float)red_ref, (float)red_tgt},
  };
  const int id_solve = dt_node_add(graph, module, "algnsurf", "solve",
      DT_LOCAL_SIZE_X, DT_LOCAL_SIZE_Y, 1, sizeof(spc), (const int *)&spc, 4,
      "ref",   "read",  "ssbo", "f32", dt_no_roi,
      "tgt",   "read",  "ssbo", "f32", dt_no_roi,
      "match", "read",  "ssbo", "u32", dt_no_roi,
      "h",     "write", "ssbo", "u32", &h_roi);
  const int id_align = dt_node_add(graph, module, "algnsurf", "align",
      roi_out.wd, roi_out.ht, 1, 0, 0, 3,
      "input2", "read",  "rgba", "*",   dt_no_roi,
      "output", "write", "rgba", "f16", &roi_out,
      "h",      "read",  "ssbo", "u32", dt_no_roi);

  CONN(dt_node_connect(graph, id_ref, 2, id_match, 0));
  CONN(dt_node_connect(graph, id_tgt, 2, id_match, 1));
  CONN(dt_node_connect(graph, id_ref, 2, id_solve, 0));
  CONN(dt_node_connect(graph, id_tgt, 2, id_solve, 1));
  CONN(dt_node_connect(graph, id_match, 2, id_solve, 2));
  CONN(dt_node_connect(graph, id_solve, 3, id_align, 2));
  dt_connector_copy(graph, module, 1, id_align, 0);
  dt_connector_copy(graph, module, 2, id_align, 1);
  graph->node[id_align].connector[0].flags |= s_conn_smooth; // sample the target with bilinear filtering
}
