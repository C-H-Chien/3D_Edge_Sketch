import os
import numpy as np
import argparse
import open3d as o3d
import logging
import sys

from eval_util import (
    set_random_seeds,
    compute_chamfer_distance,
    f_score,
    compute_precision_recall_IOU
)

def txt_to_numpy_array(file_path, delimiter=None, dtype=None):
    try:
        # Use np.loadtxt to read the data from the text file
        array = np.loadtxt(file_path, delimiter=delimiter, dtype=dtype)
        return array
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

def update_totals_and_metrics(metrics, totals, results, edge_type):
    correct_gt, num_gt, correct_pred, num_pred, acc, comp = results
    metrics[f"comp_{edge_type}"].append(comp)
    metrics[f"acc_{edge_type}"].append(acc)
    for i, threshold in enumerate(["5", "10", "20"]):
        totals[f"thre{threshold}_correct_gt_total"] += correct_gt[i]
        totals[f"thre{threshold}_correct_pred_total"] += correct_pred[i]
    totals["num_gt_total"] += num_gt
    totals["num_pred_total"] += num_pred


def finalize_metrics(metrics):
    for key, value in metrics.items():
        value = np.array(value)
        value[np.isnan(value)] = 0
        metrics[key] = round(np.mean(value), 4)
    return metrics


def print_metrics(metrics, totals, edge_type):
    print(f"{edge_type.capitalize()}:")
    print(f"  Completeness: {metrics[f'comp_{edge_type}']}")
    print(f"  Accuracy: {metrics[f'acc_{edge_type}']}")


def eval_3D_edges(obj_name, base_dir, output_name, GT_dir, metrics, totals):
    #> get the 3D edge sketch curve points
    print(f"Processing: {obj_name}")
    out_cruve_points_path = os.path.join( base_dir, output_name, "All_3D_edges_ABC-NEF_00004605.txt" )
    print( out_cruve_points_path )
    if not os.path.exists( out_cruve_points_path ):
        print(f"Invalid 3D edge sketch result at {obj_name}")
        return

    edge_sketch_curve_points = txt_to_numpy_array( out_cruve_points_path )
    if len(edge_sketch_curve_points) == 0:
        print(f"Invalid 3D edge sketch at {obj_name}")
        return

    print(edge_sketch_curve_points)

    #> get the ground-truth 3D curve points
    gt_cruve_points_path = os.path.join( base_dir, GT_dir, f'{obj_name}.txt' )
    print( gt_cruve_points_path )
    if not os.path.exists( gt_cruve_points_path ):
        print(f"Invalid GT curve points at {gt_cruve_points_path}")
        return
    
    gt_points = txt_to_numpy_array( gt_cruve_points_path )
    if len(gt_points) == 0:
        print(f"Invalid GT curve points at {obj_name}")
        return

    chamfer_dist, acc, comp = compute_chamfer_distance(edge_sketch_curve_points, gt_points)
    print(
        f"  Chamfer Distance: {chamfer_dist:.4f}, Accuracy: {acc:.4f}, Completeness: {comp:.4f}"
    )
    metrics["chamfer"].append(chamfer_dist)
    metrics["acc"].append(acc)
    metrics["comp"].append(comp)

    # #> Write GT edge points to a ply file
    # gt_edge_pcd = o3d.geometry.PointCloud()
    # gt_edge_pcd.points = o3d.utility.Vector3dVector(gt_points)

    # gt_edge_ply_file_path = os.path.join(base_dir, obj_name, exp_name, "results", "gt_edge_points.ply")
    # try:
    #     o3d.io.write_point_cloud(gt_edge_ply_file_path, gt_edge_pcd, write_ascii=True)
    #     logging.info(f"Saved {gt_edge_ply_file_path} for edge points visualization.")
    # except IOError as e:
    #     logging.error(f"Failed to save {gt_edge_ply_file_path}: {e}")

    metrics = compute_precision_recall_IOU(
        edge_sketch_curve_points,
        gt_points,
        metrics,
        thresh_list=[0.005, 0.01, 0.02],
        edge_type="all",
    )


def main(base_dir, GT_dir, exp_name):
    eval_obj_name = "00004605"
    set_random_seeds()
    metrics = {
        "chamfer": [],
        "acc": [],
        "comp": [],
        "comp_curve": [],
        "comp_line": [],
        "acc_curve": [],
        "acc_line": [],
        "precision_0.01": [],
        "recall_0.01": [],
        "fscore_0.01": [],
        "IOU_0.01": [],
        "precision_0.02": [],
        "recall_0.02": [],
        "fscore_0.02": [],
        "IOU_0.02": [],
        "precision_0.005": [],
        "recall_0.005": [],
        "fscore_0.005": [],
        "IOU_0.005": [],
    }

    totals = {
        "curve": {
            "thre5_correct_gt_total": 0,
            "thre10_correct_gt_total": 0,
            "thre20_correct_gt_total": 0,
            "thre5_correct_pred_total": 0,
            "thre10_correct_pred_total": 0,
            "thre20_correct_pred_total": 0,
            "num_gt_total": 0,
            "num_pred_total": 0,
        },
        "line": {
            "thre5_correct_gt_total": 0,
            "thre10_correct_gt_total": 0,
            "thre20_correct_gt_total": 0,
            "thre5_correct_pred_total": 0,
            "thre10_correct_pred_total": 0,
            "thre20_correct_pred_total": 0,
            "num_gt_total": 0,
            "num_pred_total": 0,
        },
    }

    eval_3D_edges(eval_obj_name, base_dir, exp_name, GT_dir, metrics, totals)

    metrics = finalize_metrics(metrics)

    print("Summary:")
    print(f"  Accuracy: {metrics['acc']:.4f}")
    print(f"  Completeness: {metrics['comp']:.4f}")
    print(f"  Recall @ 5 mm: {metrics['recall_0.005']:.4f}")
    print(f"  Recall @ 10 mm: {metrics['recall_0.01']:.4f}")
    print(f"  Recall @ 20 mm: {metrics['recall_0.02']:.4f}")
    print(f"  Precision @ 5 mm: {metrics['precision_0.005']:.4f}")
    print(f"  Precision @ 10 mm: {metrics['precision_0.01']:.4f}")
    print(f"  Precision @ 20 mm: {metrics['precision_0.02']:.4f}")
    print(f"  F-Score @ 5 mm: {metrics['fscore_0.005']:.4f}")
    print(f"  F-Score @ 10 mm: {metrics['fscore_0.01']:.4f}")
    print(f"  F-Score @ 20 mm: {metrics['fscore_0.02']:.4f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process CAD data and compute metrics."
    )
    parser.add_argument(
        "--base_dir",
        type=str,
        default="/gpfs/data/bkimia/cchien3/3D_Edge_Sketch/",
        help="Base directory for 3D Edge Sketch",
    )
    parser.add_argument(
        "--GT_dir",
        type=str,
        default="/gpfs/data/bkimia/Datasets/ABC-NEF/gt_curve_points/",
        help="Directory for the GT 3D curve points",
    )
    parser.add_argument("--output_name", type=str, default="outputs", help="output folder name")

    args = parser.parse_args()
    main(args.base_dir, args.GT_dir, args.output_name)
