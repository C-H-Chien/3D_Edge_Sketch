import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import colorsys
import random
import json

def visualize_gt(all_gt_points, name, save_fig=False, show_fig=True):
    ax = plt.figure(dpi=120).add_subplot(projection='3d')
    x = [k[0] for k in all_gt_points]
    y = [k[1] for k in all_gt_points]
    z = [k[2] for k in all_gt_points]
    ax.scatter(x, y, z, c='g', marker='o', s=0.5, linewidth=1, alpha=1, cmap='spectral')

    # ax.axis('auto')
    plt.axis('off')
    # plt.xlabel("X axis")
    # plt.ylabel("Y axis")

    with open(os.path.join(vis_gt_dir, name + ".txt"), "w") as file:
        for k in all_gt_points:
            file.write('{}'.format(k[0]))
            file.write('\t{}'.format(k[1]))
            file.write('\t{}\n'.format(k[2]))

    ax.view_init(azim=60, elev=60)
    range_size = [0, 1]
    ax.set_zlim3d(range_size[0], range_size[1])
    plt.axis([range_size[0], range_size[1], range_size[0], range_size[1]])
    if save_fig:
        plt.savefig(os.path.join(vis_gt_dir, name + ".png"), bbox_inches='tight')
    if show_fig:
        plt.show()

def get_gt_points(name, base_dir):
    objs_dir = os.path.join(base_dir, "obj")
    obj_names = os.listdir(objs_dir)
    obj_names.sort()
    index_obj_names = {}
    for obj_name in obj_names:
        index_obj_names[obj_name[:8]] = obj_name

    json_feats_path = os.path.join(base_dir, "chunk_0000_feats.json")
    with open(json_feats_path, 'r') as f:
        json_data_feats = json.load(f)
    json_stats_path = os.path.join(base_dir, "chunk_0000_stats.json")
    with open(json_stats_path, 'r') as f:
        json_data_stats = json.load(f)

    # get the normalize scale to help align the nerf points and gt points
    [x_min, y_min, z_min, x_max, y_max, z_max, x_range, y_range, z_range] = json_data_stats[name]["bbox"]
    scale = 1 / max(x_range, y_range, z_range)
    # print("normalize scale:", scale)
    poi_center = np.array([((x_min + x_max) / 2), ((y_min + y_max) / 2), ((z_min + z_max) / 2)]) * scale
    # print("poi:", poi_center)
    set_location = [0.5, 0.5, 0.5] - poi_center  # based on the rendering settings

    obj_path = os.path.join(objs_dir, index_obj_names[name])
    with open(obj_path, encoding='utf-8') as file:
        data = file.readlines()
    vertices_obj = [each.split(' ') for each in data if each.split(' ')[0] == 'v']
    vertices_xyz = [[float(v[1]), float(v[2]), float(v[3].replace('\n', ''))] for v in vertices_obj]

    edge_pts = []
    edge_pts_raw = []
    for each_curve in json_data_feats[name]:
        #> Take all curves, regardless whether any of them is 'sharp' or not
        # if each_curve['sharp']:
        each_edge_pts = [vertices_xyz[i] for i in each_curve['vert_indices']]
        edge_pts_raw += each_edge_pts

        gt_sampling = []
        each_edge_pts = np.array(each_edge_pts)
        for index in range(len(each_edge_pts) - 1):
            next = each_edge_pts[index + 1]
            current = each_edge_pts[index]
            num = int(np.linalg.norm(next - current) // 0.01)
            linspace = np.linspace(0, 1, num)
            gt_sampling.append(linspace[:, None] * current + (1 - linspace)[:, None] * next)
        each_edge_pts = np.concatenate(gt_sampling).tolist()
        edge_pts += each_edge_pts


    edge_pts_raw = np.array(edge_pts_raw) * scale + set_location
    edge_pts = np.array(edge_pts) * scale + set_location

    return edge_pts_raw.astype(np.float32), edge_pts.astype(np.float32)

#> Define the root paths and base directory
root_dir = "/gpfs/data/bkimia/cchien3/Qiwu_Edge_Grouping/3D_Edge_Sketch_and_Grouping/"
base_dir = root_dir + "ABC_NEF_obj"
objs_dir = os.path.join(base_dir, "obj")
obj_names = os.listdir(objs_dir)
obj_names.sort()
sorted_obj_names = np.array(obj_names)
for i, obj_name in enumerate(obj_names):
    sorted_obj_names[i] = obj_name[:8]

#> Define the path for saving figures
vis_gt_dir = root_dir + "evaluations/gt_curve_points"
os.makedirs(vis_gt_dir, exist_ok=True)
desired_obj = "00004605"
do_only_on_desired_obj = True

if do_only_on_desired_obj == True:
    #> find the desired object in the sorted_obj_names
    for i, name in enumerate(sorted_obj_names):
        if name == desired_obj:
            break

    print("-" * 50)
    print("processing:", i, ", name:", name)

    gt_points_raw, gt_points = get_gt_points(name, base_dir)
    visualize_gt(gt_points, name, save_fig=True, show_fig=True)
else:
    for i, name in enumerate(sorted_obj_names):
        print("-" * 50)
        print("processing:", i, ", name:", name)

        gt_points_raw, gt_points = get_gt_points(name, base_dir)
        visualize_gt(gt_points, name, save_fig=True, show_fig=False)