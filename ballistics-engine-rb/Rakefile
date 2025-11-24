require "bundler/gem_tasks"
require "rb_sys/extensiontask"

# Load gemspec
spec = Gem::Specification.load("ballistics-engine.gemspec")

# Build the Rust extension with rb-sys
RbSys::ExtensionTask.new("ballistics_engine_rb", spec) do |ext|
  ext.lib_dir = "lib"
end

task default: [:compile]
