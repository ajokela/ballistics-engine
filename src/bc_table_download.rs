//! BC5D Table Auto-Download Module
//!
//! This module provides automatic downloading and caching of BC5D correction tables
//! from a remote server. Tables are cached locally to avoid repeated downloads.
//!
//! # Example
//!
//! ```no_run
//! use ballistics_engine::bc_table_download::Bc5dDownloader;
//!
//! let downloader = Bc5dDownloader::new(
//!     "https://ballistics.tools/downloads/bc5d",
//!     false
//! ).unwrap();
//!
//! // Download table for .308 caliber
//! let table_path = downloader.ensure_table(0.308).unwrap();
//! ```

use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::PathBuf;
use std::time::Duration;

/// Default URL for BC5D table downloads
pub const DEFAULT_BC5D_URL: &str = "https://ballistics.tools/downloads/bc5d";

/// Download timeout in seconds
const DOWNLOAD_TIMEOUT_SECS: u64 = 60;

/// Manifest file name
const MANIFEST_FILE: &str = "manifest.json";

/// Error type for BC5D table download operations
#[derive(Debug)]
pub enum Bc5dDownloadError {
    /// Network error during download
    NetworkError(String),
    /// Request timed out
    Timeout,
    /// IO error reading/writing files
    IoError(std::io::Error),
    /// CRC32 checksum mismatch after download
    ChecksumMismatch { expected: String, actual: String },
    /// Requested caliber not available
    CaliberNotAvailable { requested: f64, available: Vec<f64> },
    /// Failed to parse manifest
    ManifestParseError(String),
    /// Cache directory could not be created
    CacheDirectoryError(String),
}

impl std::fmt::Display for Bc5dDownloadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Bc5dDownloadError::NetworkError(msg) => write!(f, "Network error: {}", msg),
            Bc5dDownloadError::Timeout => write!(f, "Download timed out"),
            Bc5dDownloadError::IoError(e) => write!(f, "IO error: {}", e),
            Bc5dDownloadError::ChecksumMismatch { expected, actual } => {
                write!(f, "Checksum mismatch: expected {}, got {}", expected, actual)
            }
            Bc5dDownloadError::CaliberNotAvailable { requested, available } => {
                let available_str: Vec<String> = available.iter().map(|c| format!(".{}", (c * 1000.0) as i32)).collect();
                write!(
                    f,
                    "No BC5D table available for caliber {:.3} ({:.1}mm)\nAvailable calibers: {}",
                    requested,
                    requested * 25.4,
                    available_str.join(", ")
                )
            }
            Bc5dDownloadError::ManifestParseError(msg) => write!(f, "Manifest parse error: {}", msg),
            Bc5dDownloadError::CacheDirectoryError(msg) => write!(f, "Cache directory error: {}", msg),
        }
    }
}

impl std::error::Error for Bc5dDownloadError {}

impl From<std::io::Error> for Bc5dDownloadError {
    fn from(e: std::io::Error) -> Self {
        Bc5dDownloadError::IoError(e)
    }
}

/// Manifest entry for a BC5D table
#[derive(Debug, Clone)]
pub struct TableEntry {
    /// Filename
    pub file: String,
    /// File size in bytes
    pub size: u64,
    /// CRC32 checksum (hex string)
    pub crc32: String,
}

/// Manifest describing available BC5D tables
#[derive(Debug, Clone)]
pub struct Bc5dManifest {
    /// Manifest version
    pub version: String,
    /// When the manifest was generated
    pub generated: String,
    /// Map of caliber (as string like "308") to table entry
    pub tables: std::collections::HashMap<String, TableEntry>,
}

/// BC5D table downloader with caching
pub struct Bc5dDownloader {
    /// Base URL for downloads
    base_url: String,
    /// Local cache directory
    cache_dir: PathBuf,
    /// Force re-download even if cached
    force_refresh: bool,
    /// Cached manifest
    manifest: Option<Bc5dManifest>,
}

impl Bc5dDownloader {
    /// Create a new downloader
    ///
    /// # Arguments
    /// * `base_url` - Base URL for BC5D table downloads
    /// * `force_refresh` - If true, always re-download even if cached
    ///
    /// # Returns
    /// A new downloader, or an error if the cache directory cannot be created
    pub fn new(base_url: &str, force_refresh: bool) -> Result<Self, Bc5dDownloadError> {
        let cache_dir = get_cache_directory()?;

        // Create cache directory if it doesn't exist
        if !cache_dir.exists() {
            fs::create_dir_all(&cache_dir).map_err(|e| {
                Bc5dDownloadError::CacheDirectoryError(format!(
                    "Failed to create cache directory {}: {}",
                    cache_dir.display(),
                    e
                ))
            })?;
        }

        Ok(Bc5dDownloader {
            base_url: base_url.trim_end_matches('/').to_string(),
            cache_dir,
            force_refresh,
            manifest: None,
        })
    }

    /// Ensure a BC5D table is available for the given caliber
    ///
    /// Downloads the table if not cached or if force_refresh is true.
    ///
    /// # Arguments
    /// * `caliber` - Bullet caliber in inches (e.g., 0.308)
    ///
    /// # Returns
    /// Path to the cached table file
    pub fn ensure_table(&mut self, caliber: f64) -> Result<PathBuf, Bc5dDownloadError> {
        // Load manifest if not already loaded
        if self.manifest.is_none() {
            self.manifest = Some(self.fetch_manifest()?);
        }
        let manifest = self.manifest.as_ref().unwrap();

        // Convert caliber to key (e.g., 0.308 -> "308")
        let caliber_key = format!("{}", (caliber * 1000.0).round() as i32);

        // Check if caliber is available
        let entry = manifest.tables.get(&caliber_key).ok_or_else(|| {
            Bc5dDownloadError::CaliberNotAvailable {
                requested: caliber,
                available: self.available_calibers_from_manifest(manifest),
            }
        })?;

        // Check if already cached and valid
        let cached_path = self.cache_dir.join(&entry.file);
        if !self.force_refresh && cached_path.exists() {
            // Verify checksum
            if let Ok(actual_crc) = calculate_file_crc32(&cached_path) {
                if actual_crc == entry.crc32 {
                    return Ok(cached_path);
                }
                // Checksum mismatch - re-download
                eprintln!("Warning: Cached table checksum mismatch, re-downloading...");
            }
        }

        // Download the table
        self.download_table(&entry.file, &cached_path, &entry.crc32)?;

        Ok(cached_path)
    }

    /// Get list of available calibers from the server
    pub fn available_calibers(&mut self) -> Result<Vec<f64>, Bc5dDownloadError> {
        if self.manifest.is_none() {
            self.manifest = Some(self.fetch_manifest()?);
        }
        Ok(self.available_calibers_from_manifest(self.manifest.as_ref().unwrap()))
    }

    /// Get the cache directory path
    pub fn cache_dir(&self) -> &PathBuf {
        &self.cache_dir
    }

    /// Check if a table is cached for the given caliber
    pub fn is_cached(&self, caliber: f64) -> bool {
        let caliber_key = (caliber * 1000.0).round() as i32;
        let filename = format!("bc5d_{}.bin", caliber_key);
        self.cache_dir.join(&filename).exists()
    }

    /// Get available calibers from a manifest
    fn available_calibers_from_manifest(&self, manifest: &Bc5dManifest) -> Vec<f64> {
        let mut calibers: Vec<f64> = manifest
            .tables
            .keys()
            .filter_map(|k| k.parse::<i32>().ok())
            .map(|k| k as f64 / 1000.0)
            .collect();
        calibers.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        calibers
    }

    /// Fetch the manifest from the server
    #[cfg(feature = "online")]
    fn fetch_manifest(&self) -> Result<Bc5dManifest, Bc5dDownloadError> {
        let url = format!("{}/{}", self.base_url, MANIFEST_FILE);

        let response = ureq::get(&url)
            .timeout(Duration::from_secs(DOWNLOAD_TIMEOUT_SECS))
            .call()
            .map_err(|e| match e {
                ureq::Error::Transport(t) if t.kind() == ureq::ErrorKind::Io => {
                    Bc5dDownloadError::NetworkError(format!("Connection failed: {}", t))
                }
                _ => Bc5dDownloadError::NetworkError(format!("{}", e)),
            })?;

        let json: serde_json::Value = response.into_json().map_err(|e| {
            Bc5dDownloadError::ManifestParseError(format!("Failed to parse JSON: {}", e))
        })?;

        // Parse manifest
        let version = json["version"]
            .as_str()
            .unwrap_or("unknown")
            .to_string();
        let generated = json["generated"]
            .as_str()
            .unwrap_or("unknown")
            .to_string();

        let tables_obj = json["tables"]
            .as_object()
            .ok_or_else(|| Bc5dDownloadError::ManifestParseError("Missing 'tables' field".to_string()))?;

        let mut tables = std::collections::HashMap::new();
        for (caliber, entry) in tables_obj {
            let file = entry["file"]
                .as_str()
                .ok_or_else(|| Bc5dDownloadError::ManifestParseError(format!("Missing 'file' for caliber {}", caliber)))?
                .to_string();
            let size = entry["size"]
                .as_u64()
                .ok_or_else(|| Bc5dDownloadError::ManifestParseError(format!("Missing 'size' for caliber {}", caliber)))?;
            let crc32 = entry["crc32"]
                .as_str()
                .ok_or_else(|| Bc5dDownloadError::ManifestParseError(format!("Missing 'crc32' for caliber {}", caliber)))?
                .to_string();

            tables.insert(caliber.clone(), TableEntry { file, size, crc32 });
        }

        Ok(Bc5dManifest {
            version,
            generated,
            tables,
        })
    }

    /// Stub for non-online builds
    #[cfg(not(feature = "online"))]
    fn fetch_manifest(&self) -> Result<Bc5dManifest, Bc5dDownloadError> {
        Err(Bc5dDownloadError::NetworkError(
            "Online features not enabled. Build with --features online".to_string(),
        ))
    }

    /// Download a table file
    #[cfg(feature = "online")]
    fn download_table(&self, filename: &str, dest_path: &PathBuf, expected_crc: &str) -> Result<(), Bc5dDownloadError> {
        let url = format!("{}/{}", self.base_url, filename);

        eprintln!("Downloading BC5D table: {}...", filename);

        let response = ureq::get(&url)
            .timeout(Duration::from_secs(DOWNLOAD_TIMEOUT_SECS))
            .call()
            .map_err(|e| match e {
                ureq::Error::Transport(t) if t.kind() == ureq::ErrorKind::Io => {
                    Bc5dDownloadError::NetworkError(format!("Connection failed: {}", t))
                }
                _ => Bc5dDownloadError::NetworkError(format!("{}", e)),
            })?;

        // Read response body
        let mut data = Vec::new();
        response.into_reader().read_to_end(&mut data).map_err(|e| {
            Bc5dDownloadError::NetworkError(format!("Failed to read response: {}", e))
        })?;

        // Verify checksum
        let actual_crc = calculate_crc32(&data);
        if actual_crc != expected_crc {
            return Err(Bc5dDownloadError::ChecksumMismatch {
                expected: expected_crc.to_string(),
                actual: actual_crc,
            });
        }

        // Write to file
        let mut file = File::create(dest_path)?;
        file.write_all(&data)?;

        eprintln!("Downloaded {} ({} bytes)", filename, data.len());

        Ok(())
    }

    /// Stub for non-online builds
    #[cfg(not(feature = "online"))]
    fn download_table(&self, _filename: &str, _dest_path: &PathBuf, _expected_crc: &str) -> Result<(), Bc5dDownloadError> {
        Err(Bc5dDownloadError::NetworkError(
            "Online features not enabled. Build with --features online".to_string(),
        ))
    }
}

/// Get the platform-specific cache directory for BC5D tables
pub fn get_cache_directory() -> Result<PathBuf, Bc5dDownloadError> {
    // Try to use the dirs crate for platform-specific directories
    if let Some(cache_dir) = dirs::cache_dir() {
        return Ok(cache_dir.join("ballistics-engine").join("bc5d"));
    }

    // Fallback to home directory
    if let Some(home) = dirs::home_dir() {
        #[cfg(target_os = "macos")]
        return Ok(home.join("Library").join("Caches").join("ballistics-engine").join("bc5d"));

        #[cfg(target_os = "windows")]
        return Ok(home.join("AppData").join("Local").join("ballistics-engine").join("cache").join("bc5d"));

        #[cfg(not(any(target_os = "macos", target_os = "windows")))]
        return Ok(home.join(".cache").join("ballistics-engine").join("bc5d"));
    }

    Err(Bc5dDownloadError::CacheDirectoryError(
        "Could not determine cache directory".to_string(),
    ))
}

/// Calculate CRC32 of a byte slice
fn calculate_crc32(data: &[u8]) -> String {
    const TABLE: [u32; 256] = make_crc32_table();
    let mut crc = 0xFFFFFFFFu32;
    for &byte in data {
        let idx = ((crc ^ byte as u32) & 0xFF) as usize;
        crc = (crc >> 8) ^ TABLE[idx];
    }
    format!("{:08x}", !crc)
}

/// Calculate CRC32 of a file
fn calculate_file_crc32(path: &PathBuf) -> Result<String, std::io::Error> {
    let mut file = File::open(path)?;
    let mut data = Vec::new();
    file.read_to_end(&mut data)?;
    Ok(calculate_crc32(&data))
}

/// Generate CRC32 lookup table (IEEE polynomial)
const fn make_crc32_table() -> [u32; 256] {
    const POLY: u32 = 0xEDB88320;
    let mut table = [0u32; 256];
    let mut i = 0;
    while i < 256 {
        let mut crc = i as u32;
        let mut j = 0;
        while j < 8 {
            if crc & 1 != 0 {
                crc = (crc >> 1) ^ POLY;
            } else {
                crc >>= 1;
            }
            j += 1;
        }
        table[i] = crc;
        i += 1;
    }
    table
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crc32_calculation() {
        // Test with known value
        let data = b"123456789";
        let crc = calculate_crc32(data);
        assert_eq!(crc, "cbf43926");
    }

    #[test]
    fn test_cache_directory() {
        let cache_dir = get_cache_directory();
        assert!(cache_dir.is_ok());
        let path = cache_dir.unwrap();
        assert!(path.to_string_lossy().contains("bc5d"));
    }

    #[test]
    fn test_caliber_key_conversion() {
        // Test caliber to key conversion
        let caliber: f64 = 0.308;
        let key = format!("{}", (caliber * 1000.0).round() as i32);
        assert_eq!(key, "308");

        let caliber: f64 = 0.224;
        let key = format!("{}", (caliber * 1000.0).round() as i32);
        assert_eq!(key, "224");
    }

    #[test]
    fn test_error_display() {
        let err = Bc5dDownloadError::CaliberNotAvailable {
            requested: 0.375,
            available: vec![0.224, 0.308, 0.338],
        };
        let msg = format!("{}", err);
        assert!(msg.contains("0.375"));
        assert!(msg.contains("9.5mm"));
        assert!(msg.contains(".224"));
        assert!(msg.contains(".308"));
        assert!(msg.contains(".338"));
    }
}
