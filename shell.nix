{
  system ? builtins.currentSystem,
}:

let
  dev = import ./dependencies.nix { inherit system; };
in

dev.pkgs.mkShell {

  buildInputs =
    with dev.pkgs;
    [
      dev.rustToolchain
    ]
    ++ dev.devTools
    ++ dev.platformDeps;

  shellHook = ''
    export TMPDIR="/tmp/nix-shell-$$"
    mkdir -p "$TMPDIR"

    echo "Platform: ${dev.pkgs.stdenv.system}"
    echo "Rust: $(rustc --version)"
    echo "Cargo: $(cargo --version)"
    echo ""
    echo "All tools pinned to nixpkgs 25.05"
  '';

  # RUST_BACKTRACE = "full";
  # RUST_LOG = "debug";
}
