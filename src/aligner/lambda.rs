use std::any::Any;
use std::os::raw::{c_int, c_void};
use std::panic::{self, AssertUnwindSafe};
use std::sync::Mutex;

pub(crate) struct LambdaMatcherContext<'a, F>
where
    F: Fn(usize, usize) -> bool + Sync,
{
    matcher: &'a F,
    panic_payload: Mutex<Option<Box<dyn Any + Send + 'static>>>,
}

impl<'a, F> LambdaMatcherContext<'a, F>
where
    F: Fn(usize, usize) -> bool + Sync,
{
    pub(crate) fn new(matcher: &'a F) -> Self {
        Self {
            matcher,
            panic_payload: Mutex::new(None),
        }
    }

    pub(crate) fn record_panic(&self, payload: Box<dyn Any + Send + 'static>) {
        let Ok(mut stored_payload) = self.panic_payload.lock() else {
            return;
        };

        if stored_payload.is_none() {
            *stored_payload = Some(payload);
        }
    }

    pub(crate) fn take_panic(&self) -> Option<Box<dyn Any + Send + 'static>> {
        self.panic_payload
            .lock()
            .ok()
            .and_then(|mut stored_payload| stored_payload.take())
    }
}

pub(crate) unsafe extern "C" fn lambda_match_trampoline<F>(
    pattern_pos: c_int,
    text_pos: c_int,
    arguments: *mut c_void,
) -> c_int
where
    F: Fn(usize, usize) -> bool + Sync,
{
    if arguments.is_null() {
        return 0;
    }

    let Ok(pattern_pos) = usize::try_from(pattern_pos) else {
        return 0;
    };
    let Ok(text_pos) = usize::try_from(text_pos) else {
        return 0;
    };

    let context = unsafe { &*(arguments as *const LambdaMatcherContext<'_, F>) };
    match panic::catch_unwind(AssertUnwindSafe(|| {
        (context.matcher)(pattern_pos, text_pos)
    })) {
        Ok(true) => 1,
        Ok(false) => 0,
        Err(payload) => {
            context.record_panic(payload);
            0
        }
    }
}
