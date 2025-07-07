use core::fmt;
use std::io;

pub type Result<T> = core::result::Result<T, Error>;

#[derive(Debug)]
pub enum Error {
    Custom(String),

    // Core
    MeshGeneration(String),

    // Morph
    Deformation(String),

    // Io
    FileSystem(io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::Custom(msg) => write!(f, "🔴 Error: {msg}"),
            Error::MeshGeneration(msg) => write!(f, "🔴 Mesh generation error: {msg}"),
            Error::Deformation(msg) => write!(f, "🔴 Deformation error: {msg}"),
            Error::FileSystem(err) => write!(f, "🔴 IO error: {err}"),
        }
    }
}

// Implementing the standard Error trait for better ecosystem compatibility.
impl std::error::Error for Error {
    // The `source` method provides the underlying cause of the error, if any.
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            // These variants don't wrap another error, so they have no source.
            Error::Custom(_) => None,
            Error::MeshGeneration(_) => None,
            Error::Deformation(_) => None,
            Error::FileSystem(err) => Some(err),
        }
    }
}

impl From<&str> for Error {
    fn from(value: &str) -> Self {
        Self::Custom(value.to_string())
    }
}

impl From<String> for Error {
    fn from(value: String) -> Self {
        Self::Custom(value)
    }
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Self {
        Self::FileSystem(err)
    }
}
