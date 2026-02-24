{
  system ? builtins.currentSystem,
}:

let
  rustOverlay = import (fetchTarball {
    url = "https://github.com/oxalica/rust-overlay/archive/37f8f092415b444c3bed6eda6bcbee51cee22e5d.tar.gz";
    sha256 = "0l9923csh45174qn3m1iwcqw1kjghdibn1ijvaknwkhhd0dlki9j";
  });
  nixpkgs = fetchTarball {
    url = "https://github.com/NixOS/nixpkgs/archive/refs/tags/25.05.tar.gz";
    sha256 = "1915r28xc4znrh2vf4rrjnxldw2imysz819gzhk9qlrkqanmfsxd";
  };

  pkgs = import nixpkgs {
    inherit system;
    overlays = [ rustOverlay ];
  };

  rustVersion = "1.88.0";
  rustToolchain = pkgs.rust-bin.stable.${rustVersion}.default.override {
    extensions = [ "rust-src" ];
  };

  cargoToml = builtins.fromTOML (builtins.readFile ./Cargo.toml);

  nativePlatformDeps =
    with pkgs;
    if stdenv.isDarwin then
      [
        # unknown until first time compiled
      ]
    else
      [
        pkg-config
      ];

  platformDeps =
    with pkgs;
    if stdenv.isDarwin then
      [
        # unsure as to which one(s) are necessary
        darwin.apple_sdk.frameworks.CoreFoundation
        darwin.apple_sdk.frameworks.CoreAudio
        darwin.apple_sdk.frameworks.AudioToolbox
      ]
    else
      [
        alsa-lib
        libjack2
        gnuplot
      ];

  devTools = with pkgs; [
    pkg-config
    clang
  ];
in
{
  inherit
    pkgs
    rustToolchain
    nativePlatformDeps
    platformDeps
    devTools
    cargoToml
    ;
}
