#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
"""Recompile OOC parser/tests against immutable, identical framework/support archives."""
import argparse
from pathlib import Path
import platform
import shlex
import subprocess

from run_l10_n50 import git, pinned, sha, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'build', 'framework-build', 'framework-source'):
        parser.add_argument('--'+name, type=Path, required=True)
    a = parser.parse_args()
    if platform.node() != '28h':
        raise RuntimeError('Build only on 28h')
    old_commit = git(a.framework_source, 'rev-parse', 'HEAD')
    pinned(a.framework_source, old_commit)
    if git(a.source, 'rev-parse', 'HEAD:src') != git(a.framework_source, 'rev-parse', 'HEAD:src'):
        raise RuntimeError('Framework headers or source differ')
    subtree = 'examples/ShapeOptimization/KelvinBall'
    changed = git(a.source, 'diff', '--name-only', old_commit, 'HEAD', '--', subtree).splitlines()
    if any(not name.startswith(subtree+'/reproduction/') and
           name not in (subtree+'/KelvinBallOptimization.cpp', subtree+'/README.md') for name in changed):
        raise RuntimeError('Kelvin support or checkpoint header changed')
    if a.build.exists():
        raise RuntimeError('Require a fresh build directory')
    (a.build/'bin').mkdir(parents=True)
    commands, immutable = [], {}

    def flags(path):
        values = {}
        for line in path.read_text().splitlines():
            if line.startswith(('CXX_DEFINES =', 'CXX_INCLUDES =', 'CXX_FLAGS =')):
                name, value = line.split('=', 1)
                values[name.strip()] = shlex.split(value.replace(str(a.framework_source), str(a.source)))
        return values['CXX_DEFINES']+values['CXX_INCLUDES']+values['CXX_FLAGS']

    def execute(command, cwd=None):
        commands.append(list(map(str,command)))
        subprocess.run(commands[-1], cwd=cwd, check=True)

    cwd = a.framework_build/'kelvin'
    native_flags = flags(cwd/'CMakeFiles/KelvinBall.dir/flags.make')
    obj = a.build/'KelvinBallOptimization.cpp.o'
    execute(['/usr/bin/c++', *native_flags, '-c', a.source/subtree/'KelvinBallOptimization.cpp', '-o', obj])
    link = shlex.split((cwd/'CMakeFiles/KelvinBall.dir/link.txt').read_text())
    link = [str(obj) if token.endswith('KelvinBallOptimization.cpp.o') else token for token in link]
    link[link.index('-o')+1] = str(a.build/'bin/KelvinBall')
    for token in link[1:]:
        path = (cwd/token).resolve()
        if path != obj and path.is_file() and token.endswith(('.a','.so','.o')):
            immutable[str(path)] = sha(path)
    execute(link, cwd)
    test_flags = flags(a.framework_build/'CMakeFiles/KelvinMemoryTest.dir/flags.make')
    objects = []
    tests = ['tests/unit/Rodin/Solver/MUMPSTest.cpp',
             'tests/unit/Rodin/Adaptation/SWIFTLinearSolverTest.cpp',
             'tests/unit/Rodin/Adaptation/SWIFTExampleParametersTest.cpp',
             'third-party/googletest/googletest/src/gtest-all.cc',
             'third-party/googletest/googletest/src/gtest_main.cc']
    for i, name in enumerate(tests):
        output = a.build/f'test-{i}.o'
        execute(['/usr/bin/c++', *test_flags,
                 '-I'+str(a.source/'third-party/googletest/googletest'),
                 '-c', a.source/name, '-o', output])
        objects.append(str(output))
    test_link = ['/usr/bin/c++', *objects, '-o', str(a.build/'bin/KelvinMemoryTest'),
                 *[str((cwd/token).resolve()) if token.endswith(('.a','.so','.o')) and not token.startswith('/')
                   else token for token in link[link.index('-o')+2:]], '-pthread', '-fopenmp']
    execute(test_link)
    if any(sha(Path(path)) != digest for path,digest in immutable.items()):
        raise RuntimeError('Immutable link inputs changed')
    write(a.build/'immutable-framework-provenance.json', dict(
        source=str(a.source), commit=git(a.source,'rev-parse','HEAD'),
        framework_source=str(a.framework_source), framework_commit=old_commit,
        framework_build=str(a.framework_build), framework_src_tree_identical=True,
        immutable_input_sha256=immutable, commands=commands,
        compiler_and_linker_exit_codes=0,
        binary_sha256=sha(a.build/'bin/KelvinBall')))


if __name__ == '__main__':
    main()
