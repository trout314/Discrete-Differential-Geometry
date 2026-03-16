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

    # Fixed histogram bins (determined from initial state)
    vtx_hist_init = sampler.manifold.degree_histogram(0)
    edge_hist_init = sampler.manifold.degree_histogram(1)

    # Vertex: even degrees only. Pick enough bins to cover initial range + headroom.
    nonzero_init = np.nonzero(vtx_hist_init)[0]
    max_vtx_deg = (nonzero_init[-1] + 1) if len(nonzero_init) > 0 else 20
    n_vtx_bins = max(20, max_vtx_deg // 2 + 5)  # number of even-degree bins
    vtx_bin_degrees = np.arange(1, n_vtx_bins + 1) * 2  # 2, 4, 6, ...
    vtx_max_deg = vtx_bin_degrees[-1]  # last regular bin degree

    # Edge: integer degrees. Cover initial range + headroom.
    nonzero_e_init = np.nonzero(edge_hist_init)[0]
    max_edge_deg = (nonzero_e_init[-1] + 1) if len(nonzero_e_init) > 0 else 10
    n_edge_bins = max(15, max_edge_deg + 3)
    edge_bin_degrees = np.arange(3, 3 + n_edge_bins)  # 3, 4, 5, ...
    edge_max_deg = edge_bin_degrees[-1]

    # Objective trajectory
    obj_history = []
    sweep_history = []

    # Style
    bar_color = "#4C72B0"
    line_color = "#C44E52"

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
        sweep_history.append(sw)

        elapsed = time.monotonic() - t_start
        tried = stats.total_tried
        accepted = stats.total_accepted
        accept_pct = 100.0 * accepted / tried if tried > 0 else 0.0

        # --- Vertex degree histogram (even degrees, fixed bins) ---
        ax_vtx.clear()
        # Map raw histogram into fixed bins; overflow goes to last bin
        raw_even = vtx_hist[1::2]  # indices 1,3,5,... -> degrees 2,4,6,...
        vtx_binned = np.zeros(n_vtx_bins, dtype=np.int64)
        n_copy = min(len(raw_even), n_vtx_bins - 1)
        vtx_binned[:n_copy] = raw_even[:n_copy]
        vtx_binned[-1] += raw_even[n_copy - 1:].sum() if n_copy > 0 else 0
        # Correct: last bin = overflow from n_copy onward
        vtx_binned[-1] = raw_even[n_vtx_bins - 1:].sum() if len(raw_even) >= n_vtx_bins else 0
        vtx_binned[:min(len(raw_even), n_vtx_bins - 1)] = raw_even[:n_vtx_bins - 1]
        if len(raw_even) >= n_vtx_bins:
            vtx_binned[-1] = raw_even[n_vtx_bins - 1:].sum()

        vtx_total = vtx_binned.sum()
        vtx_rel = vtx_binned / vtx_total if vtx_total > 0 else vtx_binned.astype(float)
        vtx_peak = vtx_rel.max()
        vtx_norm = vtx_rel / vtx_peak if vtx_peak > 0 else vtx_rel
        mean_vtx = mfd.mean_degree(0)
        ax_vtx.bar(vtx_bin_degrees, vtx_norm, color=bar_color, width=2.0)
        ax_vtx.axvline(mean_vtx, color="red", linewidth=1.2, linestyle="-")
        ax_vtx.set_xlabel("Vertex degree")
        ax_vtx.set_ylim(0, 1.05)
        ax_vtx.set_yticks(np.linspace(0, 1, 11))
        ax_vtx.tick_params(axis="y", labelleft=False)
        ax_vtx.set_xlim(vtx_bin_degrees[0] - 1, vtx_bin_degrees[-1] + 1)
        # Label last bin as overflow if it has content
        if vtx_binned[-1] > 0 and len(raw_even) >= n_vtx_bins:
            labels = [str(d) for d in vtx_bin_degrees]
            labels[-1] = f"{vtx_max_deg}+"
            ax_vtx.set_xticks(vtx_bin_degrees[::4])
            ax_vtx.set_xticklabels(labels[::4])
        if vtx_peak > 0:
            ax_vtx.text(0.97, 0.95,
                        f"peak={vtx_peak:.3f}\nmean={mean_vtx:.1f}\nvar={dv:.1f}",
                        transform=ax_vtx.transAxes, ha="right", va="top",
                        fontsize=9, color="#444444")

        # --- Edge degree histogram (fixed bins, starting at degree 3) ---
        ax_edge.clear()
        # Map raw histogram (0-indexed: index i = degree i+1) into bins starting at 3
        edge_binned = np.zeros(n_edge_bins, dtype=np.int64)
        raw_start = 2  # index 2 in edge_hist = degree 3
        for i in range(n_edge_bins - 1):
            raw_idx = raw_start + i
            if raw_idx < len(edge_hist):
                edge_binned[i] = edge_hist[raw_idx]
        # Overflow: everything from degree (3 + n_edge_bins - 1) onward
        overflow_start = raw_start + n_edge_bins - 1
        if overflow_start < len(edge_hist):
            edge_binned[-1] = edge_hist[overflow_start:].sum()

        edge_total = edge_binned.sum()
        edge_rel = edge_binned / edge_total if edge_total > 0 else edge_binned.astype(float)
        edge_peak = edge_rel.max()
        edge_norm = edge_rel / edge_peak if edge_peak > 0 else edge_rel
        mean_edge = mfd.mean_degree(1)
        ax_edge.bar(edge_bin_degrees, edge_norm, color=bar_color, width=1.0)
        ax_edge.axvline(mean_edge, color="red", linewidth=1.2, linestyle="-")
        ax_edge.set_xlabel("Edge degree (hinge degree)")
        ax_edge.set_ylim(0, 1.05)
        ax_edge.set_yticks(np.linspace(0, 1, 11))
        ax_edge.tick_params(axis="y", labelleft=False)
        ax_edge.set_xlim(edge_bin_degrees[0] - 0.5, edge_bin_degrees[-1] + 0.5)
        if edge_binned[-1] > 0 and len(edge_hist) >= n_edge_bins:
            labels_e = [str(d) for d in edge_bin_degrees]
            labels_e[-1] = f"{edge_max_deg}+"
            ax_edge.set_xticks(edge_bin_degrees)
            ax_edge.set_xticklabels(labels_e)
        edge_dv = mfd.degree_variance(1)
        if edge_peak > 0:
            ax_edge.text(0.97, 0.95,
                         f"peak={edge_peak:.3f}\nmean={mean_edge:.2f}\nvar={edge_dv:.1f}",
                         transform=ax_edge.transAxes, ha="right", va="top",
                         fontsize=9, color="#444444")

        # --- Objective trajectory ---
        ax_obj.clear()
        ax_obj.plot(sweep_history, obj_history, color=line_color, linewidth=1)
        ax_obj.set_xlabel("Sweeps")
        ax_obj.set_ylabel("Objective")

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

        fv_str = "(" + ", ".join(str(int(x)) for x in fv) + ")"

        p = sampler._params
        info = (
            f"Progress: {sw:.0f}/{total_sweeps} sweeps ({pct:.0f}%)\n"
            f"Elapsed: {elapsed:.1f}s   ETA: {eta}\n"
            f"\n"
            f"Facets: {nf:,}   (target {target_facets:,})\n"
            f"f-vector: {fv_str}\n"
            f"Objective: {obj:.1f}\n"
            f"Acceptance: {accept_pct:.1f}%  "
            f"({accepted:,}/{tried:,})\n"
            f"\n"
            f"Coefficients:\n"
            f"  \u03b2_volume:      {p.num_facets_coef}\n"
            f"  \u03b2_hinges:      {p.num_hinges_coef}\n"
            f"  \u03b2_hinge_var:   {p.hinge_degree_variance_coef}\n"
            f"  \u03b2_codim3_var:  {p.codim3_degree_variance_coef}\n"
            f"  hinge_target:  {p.hinge_degree_target}\n"
            f"\n"
            f"Move acceptance:\n" + "\n".join(move_lines)
        )
        ax_info.text(0.05, 0.95, info, transform=ax_info.transAxes,
                     fontsize=10, verticalalignment="top",
                     fontfamily="monospace")

        fig.subplots_adjust(left=0.06, right=0.94, wspace=0.12, hspace=0.22, top=0.95, bottom=0.06)
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
