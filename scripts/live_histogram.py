#!/usr/bin/env python3
"""Live histogram display for the manifold sampler.

Opens a matplotlib window showing degree histograms that update in
real time as the sampler runs. The window stays responsive because
sampling runs in a background thread between display updates.

Usage:
    PYTHONPATH=python python3 scripts/live_histogram.py [num_facets] [sweeps]

Examples:
    PYTHONPATH=python python3 scripts/live_histogram.py 1000 200
    PYTHONPATH=python python3 scripts/live_histogram.py 5000 500

Requires a working matplotlib GUI backend. If you get an import error:
    sudo apt install python3-pil.imagetk   # for TkAgg
    pip install PyQt6                       # for QtAgg
"""

import sys
import threading
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams


def main():
    target_facets = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    total_sweeps = int(sys.argv[2]) if len(sys.argv) > 2 else 200

    # --- Build sampler ---
    print(f"Building {target_facets}-facet manifold via ramped growth...")
    m = Manifold.standard_sphere(3)
    params = SamplerParams(
        num_facets_target=target_facets,
        num_facets_coef=0.1,
        num_hinges_coef=0.05,
        hinge_degree_variance_coef=0.0,
        codim3_degree_variance_coef=0.1,
        hinge_move_prob=0.3,
    )
    sampler = ManifoldSampler(m, params)
    step_size = max(500, target_facets // 20)
    sampler.ramped_grow(target_facets, step_size=step_size, eq_sweeps_per_step=3)
    print(f"Ready: {sampler.manifold.num_facets} facets, "
          f"{sampler.manifold.f_vector[0]} vertices")

    # --- Set up figure ---
    plt.ion()
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Manifold Sampler — Live Histograms", fontsize=14)
    ax_vtx, ax_edge, ax_obj, ax_info = axes.flat

    # Objective trajectory
    obj_history = []
    dv_history = []
    sweep_history = []

    # Create twin axis for degree variance (once, not per frame)
    ax_dv = ax_obj.twinx()

    # Style
    bar_color = "#4C72B0"
    line_color = "#C44E52"
    dv_color = "#55A868"

    sampler.reset_stats()
    sweeps_per_chunk = max(1.0, total_sweeps / 200)  # ~200 display updates
    sweeps_done = 0.0

    # Shared state for background sampling
    lock = threading.Lock()
    sampling_done = threading.Event()
    chunk_done = threading.Event()
    chunk_done.set()  # ready for first chunk

    def sample_loop():
        nonlocal sweeps_done
        while sweeps_done < total_sweeps:
            chunk_done.wait()
            chunk_done.clear()
            sw = min(sweeps_per_chunk, total_sweeps - sweeps_done)
            sampler.run(sweeps=sw)
            with lock:
                sweeps_done += sw
            chunk_done.set()
            time.sleep(0.01)  # yield to display thread
        sampling_done.set()

    worker = threading.Thread(target=sample_loop, daemon=True)
    worker.start()

    t_start = time.monotonic()

    while not sampling_done.is_set() or not chunk_done.is_set():
        chunk_done.wait(timeout=0.3)

        with lock:
            sw = sweeps_done
        mfd = sampler.manifold
        nf = mfd.num_facets
        fv = mfd.f_vector

        # Collect data
        vtx_hist = mfd.degree_histogram(0)
        edge_hist = mfd.degree_histogram(1)
        obj = sampler.current_objective
        dv = mfd.degree_variance(0)
        stats = sampler.get_stats()

        obj_history.append(obj)
        dv_history.append(dv)
        sweep_history.append(sw)

        elapsed = time.monotonic() - t_start
        tried = stats.total_tried
        accepted = stats.total_accepted
        accept_pct = 100.0 * accepted / tried if tried > 0 else 0.0

        # --- Vertex degree histogram ---
        ax_vtx.clear()
        degrees = np.arange(1, len(vtx_hist) + 1)
        vtx_total = vtx_hist.sum()
        vtx_freq = vtx_hist / vtx_total if vtx_total > 0 else vtx_hist
        ax_vtx.bar(degrees, vtx_freq, color=bar_color, width=0.8)
        ax_vtx.set_xlabel("Vertex degree")
        ax_vtx.set_ylabel("Frequency")
        ax_vtx.set_title(f"Vertex degrees  (n={int(fv[0])}, var={dv:.1f})")
        nonzero = np.nonzero(vtx_hist)[0]
        if len(nonzero) > 0:
            ax_vtx.set_xlim(nonzero[0], nonzero[-1] + 2)

        # --- Edge degree histogram ---
        ax_edge.clear()
        degrees_e = np.arange(1, len(edge_hist) + 1)
        edge_total = edge_hist.sum()
        edge_freq = edge_hist / edge_total if edge_total > 0 else edge_hist
        ax_edge.bar(degrees_e, edge_freq, color=bar_color, width=0.8)
        ax_edge.set_xlabel("Edge degree (hinge degree)")
        ax_edge.set_ylabel("Frequency")
        mean_edge = mfd.mean_degree(1)
        ax_edge.set_title(f"Edge degrees  (n={int(fv[1])}, mean={mean_edge:.2f})")
        nonzero_e = np.nonzero(edge_hist)[0]
        if len(nonzero_e) > 0:
            ax_edge.set_xlim(nonzero_e[0], nonzero_e[-1] + 2)

        # --- Objective + degree variance trajectory ---
        ax_obj.clear()
        ax_dv.clear()
        ax_obj.plot(sweep_history, obj_history, color=line_color, linewidth=1,
                    label="objective")
        ax_obj.set_xlabel("Sweeps")
        ax_obj.set_ylabel("Objective", color=line_color)
        ax_obj.tick_params(axis="y", labelcolor=line_color)
        ax_obj.set_title("Objective & vertex degree variance")

        ax_dv.plot(sweep_history, dv_history, color=dv_color, linewidth=1,
                   label="vtx deg var")
        ax_dv.set_ylabel("Vtx deg variance", color=dv_color)
        ax_dv.tick_params(axis="y", labelcolor=dv_color)

        # --- Info panel ---
        ax_info.clear()
        ax_info.axis("off")

        pct = 100.0 * sw / total_sweeps
        eta = ""
        if sw > 0:
            eta_s = elapsed * (total_sweeps - sw) / sw
            eta = f"{eta_s:.0f}s" if eta_s < 120 else f"{eta_s/60:.1f}m"

        bt = stats.bistellar_tries
        ba = stats.bistellar_accepts
        move_lines = []
        for i in range(4):
            if bt[i] > 0:
                r = 100.0 * ba[i] / bt[i]
                move_lines.append(f"  {i+1}->{4-i}: {ba[i]:,}/{bt[i]:,} ({r:.0f}%)")
        ht, ha = stats.hinge_tries, stats.hinge_accepts
        if ht > 0:
            move_lines.append(f"  hinge:  {ha:,}/{ht:,} ({100*ha/ht:.0f}%)")

        info = (
            f"Progress: {sw:.0f}/{total_sweeps} sweeps ({pct:.0f}%)\n"
            f"Elapsed: {elapsed:.1f}s   ETA: {eta}\n"
            f"\n"
            f"Facets: {nf:,}   (target {target_facets:,})\n"
            f"Objective: {obj:.1f}\n"
            f"Acceptance: {accept_pct:.1f}%  "
            f"({accepted:,}/{tried:,})\n"
            f"\n"
            f"Move acceptance:\n" + "\n".join(move_lines)
        )
        ax_info.text(0.05, 0.95, info, transform=ax_info.transAxes,
                     fontsize=10, verticalalignment="top",
                     fontfamily="monospace")

        fig.tight_layout()
        fig.canvas.draw_idle()
        fig.canvas.flush_events()
        plt.pause(0.05)

    worker.join()

    # Final update title
    fig.suptitle("Manifold Sampler — Complete", fontsize=14)
    fig.canvas.draw_idle()

    plt.ioff()
    plt.show()


if __name__ == "__main__":
    main()
