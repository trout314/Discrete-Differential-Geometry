project('discrete-differential-geometry', 'd',
  version : '0.1.0',
  license : 'GPL-2.0',
)

# Find the D runtime library directory for rpath so executables can find
# the shared phobos/druntime libraries at runtime.
dc = meson.get_compiler('d')
dc_path = run_command('which', dc.cmd_array()[0], check : true).stdout().strip()
dc_dir = run_command('dirname', dc_path, check : true).stdout().strip()

# LDC uses ../lib, DMD uses ../lib64
dc_lib_dir = dc_dir / '..' / 'lib'
dc_lib64_dir = dc_dir / '..' / 'lib64'
d_link_args = ['-L=-rpath=' + dc_lib_dir, '-L=-rpath=' + dc_lib64_dir]

# ---------------------------------------------------------------------------
# Source files
# ---------------------------------------------------------------------------

# Core library sources (no main() functions)
lib_sources = files(
  'source/algorithms.d',
  'source/factoring.d',
  'source/hashmap.d',
  'source/manifold.d',
  'source/manifold_examples.d',
  'source/manifold_moves.d',
  'source/polygons.d',
  'source/rational_extension_vector.d',
  'source/rational_number.d',
  'source/sampler.d',
  'source/simplicial_complex.d',
  'source/utility.d',
)

d_import_dirs = include_directories('source')

# ---------------------------------------------------------------------------
# Applications
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Shared library for Python/FFI bindings
# ---------------------------------------------------------------------------

shared_library('ddg_dlang',
  lib_sources + files(
    'source/ddg_capi.d',
  ),
  include_directories : d_import_dirs,
  d_args : ['-link-defaultlib-shared=false', '-O3', '-release'],
  name_prefix : '',
  install : true,
)

executable('bench_sampler',
  lib_sources + files('bench_sampler.d'),
  include_directories : d_import_dirs,
  link_args : d_link_args,
)

executable('bench_sampler_exact',
  lib_sources + files('bench_sampler.d'),
  include_directories : d_import_dirs,
  d_module_versions : ['TrackValidMoves'],
  link_args : d_link_args,
)

executable('test_valid_moves',
  lib_sources + files('test_valid_moves.d'),
  include_directories : d_import_dirs,
  d_module_versions : ['TrackValidMoves'],
  link_args : d_link_args,
)

# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------

# D's -unittest flag enables unittest blocks. The test_runner.d provides a
# trivial main(); D's runtime automatically runs all unittest blocks before it.
# Application main() functions are guarded by version(unittest){} else.
test_exe = executable('unittest',
  lib_sources + files(
    'source/test_runner.d',
  ),
  include_directories : d_import_dirs,
  d_unittest : true,
  link_args : d_link_args,
)

test('unittests', test_exe, workdir : meson.project_source_root())
