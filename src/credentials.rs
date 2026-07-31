//! CLI credential store for the online API Personal Access Token (PAT).
//!
//! Resolution precedence: the `BALLISTICS_API_TOKEN` environment variable takes precedence over
//! the credentials file at `~/.ballistics/credentials.toml` (mirroring the ToS-acceptance file
//! location). The file is written `0600` on Unix and never echoed after saving.

use std::fs;
use std::io;
use std::path::{Path, PathBuf};

const ENV_VAR: &str = "BALLISTICS_API_TOKEN";

/// Path to the credentials file (`~/.ballistics/credentials.toml`).
pub fn credentials_path() -> Option<PathBuf> {
    dirs::home_dir().map(|home| home.join(".ballistics").join("credentials.toml"))
}

fn save_token_at(path: &Path, token: &str) -> io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    fs::write(path, format!("token = \"{}\"\n", token.trim()))?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(path, fs::Permissions::from_mode(0o600))?;
    }
    Ok(())
}

fn load_token_from(path: &Path) -> Option<String> {
    let content = fs::read_to_string(path).ok()?;
    for line in content.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("token") {
            let rest = rest.trim_start().strip_prefix('=')?.trim();
            let val = rest.trim_matches('"');
            if !val.is_empty() {
                return Some(val.to_string());
            }
        }
    }
    None
}

/// Save the token to the credentials file (`0600` on Unix), creating the directory as needed.
pub fn save_token(token: &str) -> io::Result<()> {
    let path = credentials_path()
        .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "no home directory"))?;
    save_token_at(&path, token)
}

/// Resolve the token: `BALLISTICS_API_TOKEN` env var first, then the credentials file.
pub fn load_token() -> Option<String> {
    if let Ok(t) = std::env::var(ENV_VAR) {
        let t = t.trim().to_string();
        if !t.is_empty() {
            return Some(t);
        }
    }
    load_token_from(&credentials_path()?)
}

/// Remove the credentials file if present.
pub fn clear_token() -> io::Result<()> {
    if let Some(path) = credentials_path() {
        if path.exists() {
            fs::remove_file(path)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_path(tag: &str) -> PathBuf {
        std::env::temp_dir()
            .join(format!("ballistics-cred-{}-{}", tag, std::process::id()))
            .join("credentials.toml")
    }

    #[test]
    fn save_and_load_roundtrip() {
        let path = temp_path("rt");
        let _ = fs::remove_dir_all(path.parent().unwrap());
        save_token_at(&path, "bpat_abc123").unwrap();
        assert_eq!(load_token_from(&path).as_deref(), Some("bpat_abc123"));
        let _ = fs::remove_dir_all(path.parent().unwrap());
    }

    #[test]
    fn load_from_missing_file_is_none() {
        let path = temp_path("missing");
        let _ = fs::remove_dir_all(path.parent().unwrap());
        assert!(load_token_from(&path).is_none());
    }

    #[test]
    fn env_var_takes_precedence_over_file() {
        std::env::set_var(ENV_VAR, "bpat_env");
        assert_eq!(load_token().as_deref(), Some("bpat_env"));
        std::env::remove_var(ENV_VAR);
    }

    #[cfg(unix)]
    #[test]
    fn saved_file_is_0600() {
        use std::os::unix::fs::PermissionsExt;
        let path = temp_path("perm");
        let _ = fs::remove_dir_all(path.parent().unwrap());
        save_token_at(&path, "bpat_x").unwrap();
        let mode = fs::metadata(&path).unwrap().permissions().mode() & 0o777;
        assert_eq!(mode, 0o600);
        let _ = fs::remove_dir_all(path.parent().unwrap());
    }
}
