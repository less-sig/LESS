import os
import sys

Import("env")
Import("projenv")
platform = env.PioPlatform()

if "clang" in env.get("PIOENV", ""):
    # use clang for framework / libraries (env) and project sources (projenv)
    for e in [env, projenv]:
        e.Replace(
            AR="llvm-ar",
            AS="clang++",
            CC="clang",
            CXX="clang++",
            GDB="arm-none-eabi-gdb",  # to retain compatbility I guess
            OBJCOPY="llvm-objcopy",
            RANLIB="llvm-ranlib",
            SIZETOOL="llvm-size",
            SIZEPRINTCMD="$SIZETOOL -B -d $SOURCES",
        )
        # config file is from special llvm/clang build at https://github.com/ARM-software/LLVM-embedded-toolchain-for-Arm
        # also provides this clang-runtimes/armv6m_soft_nofp folder
        e.Prepend(
            CCFLAGS=["--config", "armv7em_hard_fpv4_sp_d16.cfg"]
        )  # add -v here to see more verbosity
        e.Prepend(ASFLAGS=["--config", "armv7em_hard_fpv4_sp_d16.cfg"])
        e.Prepend(
            LINKFLAGS=[
                "--config",
                "armv7em_hard_fpv4_sp_d16.cfg",
                "-Wl,--gc-sections",
            ]
        )
        # filter out nostdlib because it causes standard headers to not be found
        # also remove flags that are not working in clang / dont make sense
        filtered = [
            x
            for x in e["CCFLAGS"]
            if not x
            in [
                "-nostdlib",
                "--specs=nano.specs",
                "--specs=nosys.specs",
                "-fsingle-precision-constant",
            ]
        ]
        e.Replace(CCFLAGS=filtered)

        filtered2 = [
            x
            for x in e["LINKFLAGS"]
            if not x
            in [
                "-Wl,--gc-sections,--relax",
                "--specs=nano.specs",
                "--specs=nosys.specs",
            ]
        ]
        e.Replace(LINKFLAGS=filtered2)

        # -s option not available for llvm-ranlib
        filtered3 = [x for x in e["RANLIBFLAGS"] if x != "-s"]
        e.Replace(RANLIBFLAGS=filtered3)

        filtered4 = [
            x
            for x in e["ASFLAGS"]
            if not x
            in [
                "-nostdlib",
                "--specs=nano.specs",
                "--specs=nosys.specs",
                "-fsingle-precision-constant",
            ]
        ]
        e.Replace(ASFLAGS=filtered4)

        # add to path.. somehow PlatformIO does not do this although it's the toolchain package.
        pkg = platform.get_package("toolchain-clang")
        e.PrependENVPath(
            "PATH",
            os.path.join(pkg.path, "bin")
            if os.path.isdir(os.path.join(pkg.path, "bin"))
            else pkg.path,
        )

    # throw out GCC completely
    if "toolchain-gccarmnoneeabi" in platform.packages:
        del platform.packages["toolchain-gccarmnoneeabi"]
