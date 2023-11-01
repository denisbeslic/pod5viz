#!/usr/bin/python3

"""
Example use of pod5 to plot the signal data from pod 5 files.
"""
import argparse
from pathlib import Path
import numpy as np
import pod5 as p5
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import itertools


def plot_multiple_reads(pod5_input, num_reads=5, plot_out="all-example"):
    with p5.Reader(pod5_input) as reader:
        all_traces = []

        top5 = itertools.islice(reader.reads(), num_reads)
        for read in top5:
            print(f"Found read {read.read_id}")
            print(f"  Read has  {read.sample_count} samples")
            sample_rate = read.run_info.sample_rate
            signal = read.signal
            time = np.arange(len(signal)) / sample_rate
            trace = dict(
                type="scatter",
                name=f"{read.read_id}",
                x=time,
                y=signal,
            )
            all_traces.append(trace)

        fig = go.Figure(data=all_traces)
        fig.update_layout(xaxis=dict(rangeslider=dict(visible=True)))
        fig.update_layout(
            title=f"Signal for first X pod5 files",
            xaxis_title="Time (seconds)",
            yaxis_title="Raw signal (pA)",
            font_family="Times New Roman",
            font_color="blue",
            title_font_family="Times New Roman",
            title_font_color="black",
            legend_title_font_color="green",
        )
        fig.write_html(f"{plot_out}.html")


def plot_single_read(pod5_input, selected_read_id, plot_out="example"):
    with p5.Reader(pod5_input) as reader:
        read = next(reader.reads([selected_read_id]))
        sample_rate = read.run_info.sample_rate
        signal = read.signal
        time = np.arange(len(signal)) / sample_rate
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=time, y=signal))
        fig.update_layout(xaxis=dict(rangeslider=dict(visible=True)))
        fig.update_layout(
            title=f"Signal for pod5 file {pod5_input}",
            xaxis_title="Time (seconds)",
            yaxis_title="Raw signal (pA)",
            font_family="Times New Roman",
            font_color="blue",
            title_font_family="Times New Roman",
            title_font_color="black",
            legend_title_font_color="green",
        )

        fig.write_html(f"{plot_out}.html")

        plt.figure(figsize=(12,3))
        plt.plot(time[401:701], signal[401:701])
        plt.savefig(f"{plot_out}.jpg")
        plt.savefig(f"{plot_out}.svg")


def main():
    parser = argparse.ArgumentParser("Iterate through all read ids in an pod5 file")
    parser.add_argument("input", type=Path)
    args = parser.parse_args()
    plot_multiple_reads(args.input, num_reads=5)
    # TODO: Print additional information on graph f.e. total length, basecalled information
    selected_read_id = "76c01f6e-06eb-416f-88d9-eb6c7803afb8"
    plot_single_read(args.input, selected_read_id)


if __name__ == "__main__":
    main()
