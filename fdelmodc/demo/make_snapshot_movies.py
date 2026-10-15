#!/usr/bin/env python3
"""Create fixed-duration MP4 movies from fdelmodc SU snapshot files."""

import argparse
import math
import os
import struct
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter

TRACE_HEADER_BYTES = 240
FIGURE_DECORATION_WIDTH = 2
NS_OFFSET = 114
DT_OFFSET = 116
FLDR_OFFSET = 8
TRWF_OFFSET = 168
D1_OFFSET = 180
F1_OFFSET = 184
D2_OFFSET = 188
F2_OFFSET = 192
UNGPOW_OFFSET = 196


def positive_float(value):
    try:
        result = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number") from exc
    if not math.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError("must be a finite number greater than zero")
    return result


def snapshot_info(path):
    size = os.path.getsize(path)
    with open(path, "rb") as stream:
        header = stream.read(TRACE_HEADER_BYTES)
    if len(header) != TRACE_HEADER_BYTES:
        raise ValueError(f"{path}: file is too short to contain an SU trace header")

    ns = struct.unpack_from("<H", header, NS_OFFSET)[0]
    nx = struct.unpack_from("<h", header, TRWF_OFFSET)[0]
    if ns == 0 or nx <= 0:
        raise ValueError(f"{path}: invalid SU header (ns={ns}, traces per frame={nx})")

    trace_bytes = TRACE_HEADER_BYTES + ns * np.dtype("<f4").itemsize
    if size % trace_bytes:
        raise ValueError(f"{path}: size is not an integer number of {trace_bytes}-byte traces")
    trace_count = size // trace_bytes
    if trace_count % nx:
        raise ValueError(f"{path}: {trace_count} traces cannot be divided into frames of {nx}")

    frames = trace_count // nx
    with open(path, "rb") as stream:
        for frame in range(frames):
            stream.seek(frame * nx * trace_bytes + FLDR_OFFSET)
            fldr = struct.unpack("<i", stream.read(4))[0]
            if fldr != frame + 1:
                raise ValueError(
                    f"{path}: expected frame {frame + 1} at trace {frame * nx + 1}, "
                    f"found fldr={fldr}"
                )

    nz = ns
    dx = struct.unpack_from("<f", header, D2_OFFSET)[0]
    x0 = struct.unpack_from("<f", header, F2_OFFSET)[0]
    dz = struct.unpack_from("<f", header, D1_OFFSET)[0]
    z0 = struct.unpack_from("<f", header, F1_OFFSET)[0]
    snapshot_dt = struct.unpack_from("<H", header, DT_OFFSET)[0] / 1_000_000
    start_time = struct.unpack_from("<f", header, UNGPOW_OFFSET)[0]
    if dx <= 0 or dz <= 0 or snapshot_dt <= 0:
        raise ValueError(f"{path}: invalid spatial or time sampling in SU header")

    return {
        "trace_bytes": trace_bytes,
        "nx": nx,
        "nz": nz,
        "frames": frames,
        "x0": x0,
        "x1": x0 + (nx - 1) * dx,
        "z0": z0,
        "z1": z0 + (nz - 1) * dz,
        "dt": snapshot_dt,
        "t0": start_time,
    }


def read_frame(stream, info, frame_index):
    byte_count = info["nx"] * info["trace_bytes"]
    stream.seek(frame_index * byte_count)
    raw = stream.read(byte_count)
    if len(raw) != byte_count:
        raise ValueError(f"unexpected end of file while reading frame {frame_index + 1}")
    return np.ndarray(
        (info["nx"], info["nz"]),
        dtype="<f4",
        buffer=raw,
        offset=TRACE_HEADER_BYTES,
        strides=(info["trace_bytes"], np.dtype("<f4").itemsize),
    ).T


def write_movie(source, output, clip, duration, title):
    info = snapshot_info(source)
    fps = info["frames"] / duration
    figure_height = 4.5
    snapshot_aspect = (info["x1"] - info["x0"]) / (info["z1"] - info["z0"])
    figure_width = figure_height * snapshot_aspect + FIGURE_DECORATION_WIDTH
    writer = FFMpegWriter(
        fps=fps,
        codec="h264",
        extra_args=["-pix_fmt", "yuv420p"],
        metadata={"title": title},
    )

    fig, ax = plt.subplots(figsize=(figure_width, figure_height), constrained_layout=True)
    extent = (info["x0"], info["x1"], info["z1"], info["z0"])
    image = ax.imshow(
        np.zeros((info["nz"], info["nx"]), dtype=np.float32),
        cmap="seismic",
        vmin=-clip,
        vmax=clip,
        interpolation="nearest",
        extent=extent,
        aspect="equal",
    )
    ax.set_xlabel("Distance (m)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label="Amplitude")
    time_label = ax.text(
        0.99,
        0.02,
        "",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        color="black",
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )

    with open(source, "rb") as stream, writer.saving(fig, output, dpi=100):
        for index in range(info["frames"]):
            image.set_data(read_frame(stream, info, index))
            time_label.set_text(f"t = {info['t0'] + index * info['dt']:.3f} s")
            writer.grab_frame()

    plt.close(fig)
    print(
        f"Wrote {output}: {info['frames']} frames at {fps:.3f} fps "
        f"({duration:g} s), clip=±{clip:g}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Convert fdelmodc SU snapshots into fixed-duration MP4 movies."
    )
    parser.add_argument("--svz", help="vertical-velocity SU file (omit to skip this movie)")
    parser.add_argument("--sp", help="pressure SU file (omit to skip this movie)")
    parser.add_argument("--svz-output", default="snap_svz.mp4", help="vertical-velocity MP4")
    parser.add_argument("--sp-output", default="snap_sp.mp4", help="pressure MP4")
    parser.add_argument("--svz-clip", type=positive_float, default=1e-5, help="Vz color limit")
    parser.add_argument("--sp-clip", type=positive_float, default=100, help="pressure color limit")
    parser.add_argument("--duration", type=positive_float, default=20, help="movie length in seconds")
    args = parser.parse_args()

    movies = []
    if args.svz:
        movies.append((args.svz, args.svz_output, args.svz_clip, "Vertical velocity (Vz)"))
    if args.sp:
        movies.append((args.sp, args.sp_output, args.sp_clip, "Pressure (P)"))

    for source, _, _, _ in movies:
        if not os.path.isfile(source):
            parser.error(f"input file not found: {source}")
    if movies and not FFMpegWriter.isAvailable():
        parser.error(
            "ffmpeg is required to write MP4 files; install ffmpeg and make it available on PATH"
        )

    try:
        for source, output, clip, title in movies:
            write_movie(source, output, clip, args.duration, title)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
