#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(improper_ctypes)]
#![allow(deref_nullptr)]
#![allow(clippy::all)]
#![allow(clippy::missing_safety_doc)]
#![allow(clippy::pedantic)]
#![allow(clippy::restriction)]
#![allow(clippy::nursery)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
