Gem::Specification.new do |spec|
  spec.name          = "ballistics-engine"
  spec.version       = "0.13.4"
  spec.authors       = ["Alex Jokela"]
  spec.email         = ["email@tinycomputers.io"]

  spec.summary       = "High-performance ballistics calculations engine"
  spec.description   = "Ruby bindings for ballistics-engine - A high-performance Rust-based ballistics calculations library with 4-DOF trajectory modeling, wind deflection, atmospheric effects, and more."
  spec.homepage      = "https://github.com/ajokela/ballistics-engine"
  spec.license       = "MIT OR Apache-2.0"
  spec.required_ruby_version = ">= 2.7.0"

  spec.metadata["homepage_uri"] = spec.homepage
  spec.metadata["source_code_uri"] = spec.homepage
  spec.metadata["changelog_uri"] = "#{spec.homepage}/blob/main/CHANGELOG.md"
  spec.metadata["documentation_uri"] = "https://docs.ballistics.rs"

  # Specify which files should be added to the gem when it is released.
  spec.files = Dir["lib/**/*.rb", "README.md", "LICENSE*", "Cargo.toml", "Cargo.lock", "src/**/*.rs"]
  spec.require_paths = ["lib"]
  spec.extensions = ["Cargo.toml"]

  # Runtime dependencies
  spec.add_dependency "rb_sys", "~> 0.9"

  # Development dependencies
  spec.add_development_dependency "rake", "~> 13.0"
  spec.add_development_dependency "rake-compiler", "~> 1.2"
end
