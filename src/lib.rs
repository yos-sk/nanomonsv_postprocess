use std::io::{BufRead, BufReader};
use std::path::Path;
use std::error::Error;
use flate2::read::MultiGzDecoder;
use std::fs::File;

pub fn open_file<P: AsRef<Path>>(p: P) -> Result<Box<dyn BufRead>, Box<dyn Error>> {
    let r = File::open(p.as_ref())?;
    let ext = p.as_ref().extension();

    if ext == Some(std::ffi::OsStr::new("gz")) {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::new(gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::new(r);
        Ok(Box::new(buf_reader))
    }
}