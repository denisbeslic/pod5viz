#!/usr/bin/python3

"""
Example use of pod5 to plot the signal data from pod 5 files.
"""
from pysam import FastxFile
from pathlib import Path
import numpy as np
import os
import pod5 as p5
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import seaborn as sns
import pathlib
import plotly.graph_objects as go
import itertools
import click
import logging

logging.basicConfig(
    level=logging.INFO,
)
logger = logging.getLogger("pod5viz")



def read_fasta(path):
    with FastxFile(path) as fh:
        for entry in fh:
            yield (entry.sequence, entry.name)

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

def read_stats(reads, outdir):
    reads_l = list(reads)
    read_seq_len = [len(i[0]) for i in reads_l]
    logger.info(f" Minimum read length: {min(read_seq_len)}")
    logger.info(f" Mean read length: {sum(read_seq_len) / len(read_seq_len)}")
    logger.info(f" Maximum read length: {max(read_seq_len)}")
    ax = sns.violinplot(x=read_seq_len, saturation=0.5, log_scale=True, color=".8")
    sns.stripplot(x=read_seq_len, jitter=True, ax=ax, size=2.5)
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
    plt.savefig(os.path.join(outdir, "seaborn_plot.jpg"), dpi=300)
    return reads_l

def signal_stats(reads, pod5):
    read_names = [i[1] for i in reads]
    logger.info(f" Total number of reads: {len(read_names)}")
    with p5.Reader(pod5) as reader:
        for read_record in reader.reads(selection=read_names):
        #for read_name in read_names:
            print(read_record)
            #read = next(reader.reads([read_name]))
            sample_rate = read_record.run_info.sample_rate
            signal = read_record.signal
            time = np.arange(len(signal)) / sample_rate
            print(signal)
            exit()
            
        exit()
    
    for read_name in read_names:
        with p5.Reader(pod5) as reader:
            exit()
    #print(read_names)
    # TODO get mean signal length, possible drops
    exit()

def export_signals(reads, pod5):
    exit()


@click.command()
@click.option(
    "--pod5",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path
    ),
    help = "Path to multi pod5 file",
    required=True,
)
@click.option(
    "--fastq",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path
    ),
    help = "Path to fastq file",
)
@click.option(
    "--outdir",
    type=click.Path(
        file_okay=False, dir_okay=True, path_type=pathlib.Path
    ),
    help = "Path to output directory",
    default = "results/",
)
def main(pod5, fastq, outdir):
    logger.info("pod5viz")
    logger.info(f" pod5 = {pod5}")
    logger.info(f" fastq = {fastq}")
    logger.info(f" outdir = {outdir}")
    
    os.makedirs(outdir, exist_ok=True)  

    reads = read_fasta(fastq)
    reads = read_stats(reads, outdir)

    signal_stats(reads, pod5)

    # TODO export pictures of signals

    exit()

    plot_multiple_reads(args.input, num_reads=5)
    # TODO: Print additional information on graph f.e. total length, basecalled information
    #selected_read_id = "76c01f6e-06eb-416f-88d9-eb6c7803afb8"
    #plot_single_read(args.input, selected_read_id)


if __name__ == "__main__":
    main()
