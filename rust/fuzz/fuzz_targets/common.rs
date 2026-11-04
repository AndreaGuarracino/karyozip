/// Runs `f` with every refusal (`karyozip::bail`) turned into a `Deferred` panic and caught:
/// a refusal is the right answer to a forged input, any other panic is a finding.
pub fn refusals_ok(f: impl FnOnce() + std::panic::UnwindSafe) {
    karyozip::DEFER.with(|d| d.set(true));
    let r = std::panic::catch_unwind(f);
    karyozip::DEFER.with(|d| d.set(false));
    if let Err(e) = r {
        if e.downcast_ref::<karyozip::Deferred>().is_none() { std::panic::resume_unwind(e); }
    }
}
