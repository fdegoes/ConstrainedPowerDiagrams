# Prompt: describe the Linux machine that ran the overnight benchmarks

Copy the prompt below into a terminal session or an assistant on the Linux
machine (`inhuman`), where the overnight run was made. It collects only the
facts the paper needs for its machine description, and nothing else.

---

I need the hardware and OS description of this machine for a research paper.
Please run read-only commands and report exactly what they print. Do not change
any settings or install anything.

1. **Host name:** `hostname`
2. **Processor:** the CPU model string, from `lscpu` (field "Model name") or
   `/proc/cpuinfo`.
3. **Cores:**
   - physical cores per socket, sockets, and threads per core, from `lscpu`;
   - logical CPUs, from `nproc` and `getconf _NPROCESSORS_ONLN`.
   - If the machine is a VM or container, say so. Check `systemd-detect-virt`
     (or `lscpu | grep -i hypervisor`) and report the hypervisor if one is shown.
     Also report any CPU limit of the container, from
     `cat /sys/fs/cgroup/cpu.max` if it exists.
4. **Memory:** the total RAM, from `free -g` (the "Mem: total" column) and
   `grep MemTotal /proc/meminfo`. Report it in GB.
5. **Operating system:** the distribution and version, from
   `cat /etc/os-release` (PRETTY_NAME), and the kernel, from `uname -r`.
6. **Compiler and toolchain** used for the build, from `gcc --version`,
   `c++ --version` and `cmake --version`. Report the first line of each.
7. **Load right now:** `uptime`, so we know whether other jobs were running.

Report the results as a short list with these exact labels:

- Host:
- Processor:
- Sockets / cores per socket / threads per core:
- Logical CPUs (`nproc`):
- Virtualization or container limit (if any):
- RAM (GB):
- OS and version:
- Kernel:
- Compiler (gcc / c++ / cmake):
- Load average now:

Then give one sentence for the paper, in this form, filling in only the values
you measured:

> The timings were recorded on a Linux machine with a *<processor>* processor
> (*<logical CPUs>* logical cores) and *<RAM>* GB of RAM.

If a value cannot be read, write "not available" and say which command failed.
Do not estimate any value.

---

Paste the report back to the paper session. The values go into the machine
paragraph of `latex/results.tex`.
