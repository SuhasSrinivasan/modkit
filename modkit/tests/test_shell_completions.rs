use common::run_modkit;

mod common;

#[test]
fn test_generate_bash_completions() {
    let _out = run_modkit(&["generate-shell-completions", "bash"]).unwrap();
}

#[test]
fn test_generate_zsh_completions() {
    let _out = run_modkit(&["generate-shell-completions", "zsh"]).unwrap();
}

#[test]
fn test_generate_fish_completions() {
    let _out = run_modkit(&["generate-shell-completions", "fish"]).unwrap();
}
