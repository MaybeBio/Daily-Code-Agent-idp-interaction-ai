//! `OpenDNA` desktop UI.

mod app;
mod bridge;
mod components;
mod format;
mod ipc;

fn main() {
    console_error_panic_hook::set_once();
    leptos::mount::mount_to_body(app::App);
}
