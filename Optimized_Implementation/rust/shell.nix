with import <nixpkgs> {};
let 
  unstable_pkgs = import <nixos-unstable> {};
  rustToolchain = unstable_pkgs.rustPlatform.rustLibSrc;
in
{ pkgs ? import <nixpkgs> {} }:

stdenv.mkDerivation {
  name = "";
  src = ./.;

  buildInputs = [ 
    # just for cargo 
    clang 
    llvm 
    lldb

    #vscode-extensions.vadimcn.vscode-lldb
    unstable_pkgs.cargo
    unstable_pkgs.cargo-nextest
    unstable_pkgs.rustc
    unstable_pkgs.rust-analyzer

    # needed for IDE
    pkg-config  
    openssl
    git
  ] ++ (lib.optionals pkgs.stdenv.isLinux ([
    gdb
  ]));
  RUST_SRC_PATH = "${unstable_pkgs.rust.packages.stable.rustPlatform.rustLibSrc}";
}
