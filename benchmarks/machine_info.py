"""Machine and interpreter capture for the benchmark drivers (NV-08).

The archived timing tables carry no hardware field.  That is why the manuscript can
only say "a single Apple silicon laptop" and must not assert a CPU model or a memory
size for them: the precise specification ("Apple M5, 10 cores, 32 GB RAM, macOS 27.0")
exists in ``testing_suite/TEST_REPORT.md``, which documents the *test-suite* session,
not the timing session.

This module closes that gap for **future** runs only.  Every driver in this directory
merges :func:`machine_fields` into the rows it writes, so a re-issued table names the
platform, the architecture and the physical memory it was measured on, and a speed-up
quoted from it becomes traceable to a machine.

It changes nothing about what is measured: the values are read once, before any timing
loop, from ``platform`` and from the operating system's own memory query - no clock is
started, no repeat is added, and no published number moves.  The archived CSVs are
left exactly as they are.

Stdlib only, like the drivers that import it.
"""

from __future__ import annotations

import os
import platform
import subprocess
import sys

#: Sentinels for "this platform offers no cheap physical-memory query".
UNKNOWN_MEMORY = -1


def physical_memory_bytes() -> int:
    """Total physical memory in bytes, read from the OS, with no new dependency.

    ``SC_PHYS_PAGES * SC_PAGE_SIZE`` is the portable answer and is the *physical* page
    count, not the available one; macOS exposes it, Linux exposes it too, and both are
    also covered by the explicit fallbacks (``sysctl hw.memsize`` and ``/proc/meminfo``)
    so a platform that refuses ``sysconf`` still reports a number.
    """
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        if pages and page_size:
            return int(pages) * int(page_size)
    except (ValueError, OSError, KeyError, AttributeError):  # pragma: no cover - platform path
        pass
    if sys.platform == "darwin":
        try:
            out = subprocess.run(
                ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, check=False
            )
            if out.returncode == 0 and out.stdout.strip().isdigit():
                return int(out.stdout.strip())
        except OSError:  # pragma: no cover - platform path
            pass
    else:
        try:
            with open("/proc/meminfo", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("MemTotal:"):
                        return int(line.split()[1]) * 1024
        except OSError:  # pragma: no cover - platform path
            pass
    return UNKNOWN_MEMORY


def machine_fields() -> dict[str, object]:
    """The columns every benchmark row carries beside its timings.

    Names are stable across the Python and the R drivers (``os_platform``,
    ``os_machine``, ``cpu``, ``total_memory_bytes``) so a reader can line a ``_py``
    table up with its ``_r`` counterpart and see that both halves came from the same
    box - which is the claim the documentation makes and the archived tables cannot.
    """
    memory = physical_memory_bytes()
    return {
        "os_platform": platform.platform(),
        "os_machine": platform.machine(),
        "cpu": platform.processor() or f"{platform.cpu_count()} cores",
        "interpreter": f"{platform.python_implementation()} {platform.python_version()}",
        "total_memory_bytes": memory,
        "total_memory_gb": "" if memory == UNKNOWN_MEMORY else round(memory / 1e9, 1),
    }


def describe(fields: dict[str, object]) -> str:
    """One human-readable line, printed once by each driver before it measures."""
    return (
        f"measured on: {fields['os_platform']} | {fields['os_machine']} | "
        f"{fields['cpu']} | {fields['total_memory_gb']} GB physical RAM | "
        f"{fields['interpreter']}"
    )
