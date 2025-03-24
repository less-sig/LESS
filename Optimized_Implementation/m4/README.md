This is a template project for the STM32F4DISCOVERY board using PlatformIO,
with some built-in helper functions to have printf() output over SWO, cycle
count measurement, stack usage monitoring, an option to use the clang compiler
instead of gcc, as well as more esoteric features for performance
investigation (see Sections [Running code from RAM](#running-code-from-ram)
and [Stack pointer in CCM RAM](#stack-pointer-in-ccm-ram)).

The clock is initialized at 24 MHz and the ART accelerator is disabled to
reduce variability of execution times.

Modifications required for SWO to work heavily borrow from
https://github.com/maxgerhardt/pio-swo-demo, including its custom SWO viewer
target to monitor output directly from within Visual Studio Code. One can also
use `st-trace` from https://github.com/stlink-org/stlink instead by running
the following command:

`st-trace -c24`

A custom test runner reading data over SWO allows for the use of PlatformIO's
built-in unit testing facilities. Example usage is provided in the tests/
directory.

Modifications required to use clang as a compiler heavily borrow from
https://github.com/maxgerhardt/platformio-with-clang/. However, clang is
not supplied by PlatformIO, so it will be necessary to package your own
version. Instructions to do this are supplied in section
[Packaging clang](#packaging-clang)

# Board selection

Beyond the STM32F4DISCOVERY board, this kit is also compatible with the
STM32F411E-DISCO and STM32F412G-DISCO boards. However, note the MCUs used in
these boards do not have CCMRAM.

In order to use one of these boards in place of the STM32F4DISCOVERY, rename
the `platformio.ini` file in the root of this project to
`platformio_disco_f407vg.ini`, and rename either `platformio_disco_f411ve.ini`
(for the STM32F411E-DISCO board) or `platformio_disco_f412zg.ini` (for the
STM32F412G-DISCO board) to `platformio.ini`.

Initial support for the STM32F469I-DISCO board was added, however SWO isn't
working (no output). Thus, for the moment, this board is considered
unsupported.

# Periodic callback function

It is possible to declare a callback function that runs every 1 ms, in the
context of the SysTick IRQ handler. Its function prototype is:

```
void callback_1ms(void)
```

As it runs in an IRQ context, it should not perform any complex tasks which
take up too much CPU time.

## Default LED behavior while running unit tests

An example application of this callback is shown in `test/test.c`, where it is
used to blink one of the board's LEDs with ~1 Hz frequency. This is useful to
detect a hard fault while running tests; since the LED blinking code runs in
an IRQ context, it will keep blinking even if the user code is trapped in an
infinite loop. However, in case of a hard fault, its IRQ handler is called and
it cannot be preempted by the SysTick handler; thus, the LED stops blinking,
in the *on* or *off* state randomly depending on the exact instant the fault
occurred.

**To summarize:** blinking LEDs indicate a regular user code fault, while a
static LED (whether on or off) indicates a hard fault.

# Stack usage monitoring

Two facilities are provided to monitor maximum stack usage for the program.

The first uses a sampling-based technique: every 1 ms (assuming the default
SysTick configuration), the current value of the stack pointer is sampled, and
the maximum is updated if necessary. This version does not require
initialization and runs automatically. Whenever the user wants to read a
sample, they can use the macro `stack_usage_sample_get_max()` to get the
largest stack usage sampled until now. `stack_usage_sample_get_last()` is also
provided to read the most recent sample.

The second uses a watermarking-based technique. It must be initialized at the
start of `main()` by calling `stack_usage_scan_init(x)`, where `x` is the
maximum expected stack usage in bytes -- this technique is unable to detect
stack usages beyond the value provided. Whenever the user wants to read the
current value of the stack high watermark (i.e. the largest stack size seen
throughout the execution of the program), they can call the function
`stack_usage_scan()`.

Each technique has advantages and disadvantages. The sampling-based technique
may miss fast (< 1 ms) functions that use a lot of stack, but do not happen
to be running when the sample is taken. The watermarking-based technique may
report incorrect values if the same magic constant used by it is written to
the stack by the program.

Examples are shown for both facilities in the `main()` function of this
template.

# Another tool to analyze code size and stack usage: puncover

Although the stack usage monitoring facilities above help understand global
stack usage, by design it is impossible to drill down into individual
functions to locate the worst offenders. A tool that is able to report
per-function stack usage is puncover (https://github.com/HBehrens/puncover/).
It requires Python and can be installed with the following command:

```
pip install puncover
```

To use the tool, run the following command from the root directory of the
project (i.e. the directory containing the `platformio.ini` file):

```
puncover \
    --gcc_tools_base ~/.platformio/packages/toolchain-gccarmnoneeabi/bin/arm-none-eabi- \
    --elf-file .pio/build/disco_f407vg_gcc/firmware.elf \
    --build_dir .pio/build \
    --src_root .
```

This creates a local web server, and opens it in your system's default web
browser, presenting a graphical user interface for code size and stack usage
analysis.

**Note:** if using clang instead of gcc, replace the firmware location passed
to the `--elf-file` option with `.pio/build/disco_f407vg_clang/firmware.elf`.

For all the reporting features of this tool to work, it requires compiling the
binary with certain flags: `-fstack-usage -fcallgraph-info -g`. These are
already included in the template project (in the `platformio.ini` file) for
gcc. For clang, the flag `-fcallgraph-info` is not available and thus is not
included. Currently stack usage is not being reported when the project is
compiled with clang, but code size and static memory usage reports work.

However, note that some users have reported inconsistencies and errors in
puncover's output, so keep this in mind, and if possible confirm the reported
values by adding up stack usage by variables and arrays in the target function.

# Cycle count measurement

Macros are provided to aid with benchmarking, using the internal cycle counter
of Cortex-M cores.

To use them, one must first add `bm_decls;` to the same scope as the code to
be benchmarked. This declares variables required by the benchmarking facility.

The code to be benchmarked must be preceded with `bm_start()`, and followed by
`bm_end()`. To get the most recent cycle count, use the `bm_result()` macro.

For long-running functions (i.e. > 1 ms), interrupts such as SysTick will run
during the benchmark, adding noise to the benchmark results. For best results,
disable IRQs prior to calling `bm_start()` using `__disable_irq()`, and
reenable them after `bm_end()` using `__enable_irq()`. Obviously, this may
break code that depends on interrupts.

Example code is shown in the `main()` function of this template.

# Running code from RAM

For benchmarking reasons, it may be desirable to run certain functions from
RAM instead of Flash memory. gcc's `__attribute__((section("...")))` syntax
can be used to achieve this; the linker script defines a `.RamFunc` section
specifically for this. This is an example syntax of a function declaration:

```
__attribute__((section(".RamFunc"))) void f(void) {
    // ...
}
```

In assembly, this can be done by the replacing `.text` with `.RamFunc` in the
`.section` directive preceding the function.

However, due to certain particularities of the STM32F4, to achieve maximum
execution performance with code in RAM, bit 29 of the address of the function
pointer (whether written in C or assembly) must be cleared. A utility macro is
provided for this in `utils.h`, called `FIX_RAM_FUNC()`. It has a single
parameter, which is the name of the function whose address must be fixed. The
macro must be placed inside a function -- it is not possible to declare it in
the global scope. The effect of this macro is to create a function pointer to
the same function, with the same name but with the suffix `_fixed` appended to
it. For maximum performance, calls to the function should be made via this
function pointer rather than directly using the original function name.

The following example should make this clearer:

```
__attribute__((section(".RamFunc"))) void f(void) {
    // ...
}

int main() {
    // ...
    
    FIX_RAM_FUNC(f);

    __disable_irq();
    bm_start();

    f_fixed();

    bm_end();
    __enable_irq();

    // ...
}
```

It turns out that bit 29 is already clear in the Flash memory region, so it is
harmless to apply this macro to code placed in Flash memory. Therefore, if
desired, the macro can be applied indistinctly to any function being
benchmarked, regardless of whether it is actually placed in RAM. Note that the
macro itself does not influence whether a function is placed in Flash or RAM;
this is done by placing the function in a special section as explained above.

# Stack pointer in CCM RAM

Investigations revealed that code may run faster in the STM32F4 if the stack
pointer points to an address in CCM (core-coupled memory) RAM, which starts
at address `0x10000000` -- conventionally the stack pointer points to SRAM,
which starts at address `0x20000000`.

For this scenario, PlatformIO environments are supplied with the `_CCM_SP`
suffix, i.e. `disco_f407vg_gcc_CCM_SP` and `disco_f407vg_clang_CCM_SP`. By
selecting this environment when building and uploading the project, the stack
pointer will automatically point to CCM RAM. Environment names without the
`_CCM_SP` suffix performs the more conventional initialization pointing to
SRAM.

# Packaging clang

The following steps are required to create a clang package:

1. Download the desired release from
https://github.com/ARM-software/LLVM-embedded-toolchain-for-Arm/releases/.
2. Unpack the file somewhere in your filesystem; for demonstration we will
assume `/Users/xxx/clang`.
3. After unpacking, a new directory should be created inside
`/Users/xxx/clang` with the same name as the downloaded file, but without the
`.tar.gz` extension. For instance, considering the most recent release as of
this writing (16.0.0), this would be:

`/Users/xxx/clang/LLVMEmbeddedToolchainForArm-16.0.0-Darwin`

4. Move the contents of this directory (i.e. the `bin`, `lib` and other
files and directories) to `/Users/xxx/clang`.
5. Delete the `/Users/xxx/clang/LLVMEmbeddedToolchainForArm-16.0.0-Darwin`
directory:

`rmdir /Users/xxx/clang/LLVMEmbeddedToolchainForArm-16.0.0-Darwin`

6. Create the required `package.json` file in `/Users/xxx/clang/` with the
following contents:

```
{
  "name": "toolchain-clang",
  "version": "1.160000.230413",
  "description": "LLVM/Clang toolchain",
  "keywords": [
    "toolchain",
    "build tools",
    "compiler",
    "assembler",
    "linker",
    "preprocessor",
    "arm"
  ],
  "homepage": "https://github.com/llvm/llvm-project",
  "license": "Apache-2.0",
  "system": [
    "darwin_x86_64"
  ],
  "repository": {
    "type": "git",
    "url": "https://github.com/llvm/llvm-project.git"
  }
}
```

7. If necessary, change the `version` and `system` fields of `package.json`.
If the downloaded clang release has version `aa.b.c` and was released on
`yy/mm/dd`, use the following version string:

```
"version": "1.aa0b0c.yymmdd",
```

8. Similarly, if you use a different system, edit the `system` field. The
following appear to be the valid strings for the most common systems:

```
darwin_x86_64
darwin_arm64
linux_x86_64
windows_amd64
```

9. Pack the toolchain using the following commands (replacing the version
and system strings in the file name if necessary):

```
cd /Users/xxx/clang
tar zcvf toolchain-clang-darwin_x86_64-1.160000.230413.tar.gz *
```

10. Save the `toolchain-clang-darwin_x86_64-1.160000.230413.tar.gz` to a
directory of your choice. It is now possible to delete the files originally
inside `/Users/xxx/clang/` to save space. For the following, it will be
assumed that this was done if desired, and the packaged file was copied back
to `/Users/xxx/clang/`.

11. Edit `platformio.ini` and locate, inside the `[env:disco_f407vg_clang]`
section, the line that references `toolchain-clang`:

```
platform_packages = 
    toolchain-clang@file://...
```

12. Replace it with:

```
platform_packages = 
    toolchain-clang@file:///Users/xxx/clang/toolchain-clang-darwin_x86_64-1.160000.230413.tar.gz
