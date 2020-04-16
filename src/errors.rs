#[derive(Fail, Debug)]
pub enum Error {
    #[fail(display = "cannot find 0th root of unity")]
    ZeroethRootOfUnity,
    #[fail(display = "n should divide (q-1) for non-trivial root of unity")]
    TrivialRootOfUnity,
}
