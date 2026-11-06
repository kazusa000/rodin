#!/usr/bin/env python3
# Distributed under the Boost Software License, Version 1.0.
"""Compile recovery I/O against an immutable, already built identical framework."""
import argparse
import json
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
from run_l10_n50 import git, pinned, sha, write


def main():
    p = argparse.ArgumentParser()
    for name in ('source', 'build', 'campaign'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    if platform.node() != '28h':
        raise RuntimeError('Build only on 28h')
    old = json.loads(a.campaign.read_text())
    oldsource = Path(old['source'])
    pinned(oldsource, old['commit'])
    for stage in ('configure', 'compile'):
        if int((a.campaign.parent/stage/'exit-code.txt').read_text()) != 0:
            raise RuntimeError('Framework build was not successful')
    entry = next(s for s in old['stages'] if s['name']=='compile')
    oldbuild = Path(entry['command'][2])
    if sha(oldbuild/'bin/KelvinBall') != old['binary_sha256']:
        raise RuntimeError('Original compiled binary changed')
    changed = git(a.source, 'diff', '--name-only', old['commit'], 'HEAD', '--', 'src', 'examples').splitlines()
    optimization = 'examples/ShapeOptimization/KelvinBall/KelvinBallOptimization.cpp'
    if any(name != optimization and not name.startswith('examples/ShapeOptimization/KelvinBall/reproduction/') for name in changed):
        raise RuntimeError('Framework/support code or checkpoint header changed')
    a.build.mkdir(exist_ok=False)
    (a.build/'bin').mkdir()
    cwd = oldbuild/'kelvin'
    make = cwd/'CMakeFiles/KelvinBall.dir/flags.make'
    flags = {}
    for line in make.read_text().splitlines():
        if line.startswith(('CXX_DEFINES =', 'CXX_INCLUDES =', 'CXX_FLAGS =')):
            name, value = line.split('=',1)
            flags[name.strip()] = shlex.split(value.replace(str(oldsource),str(a.source)))
    obj = a.build/'KelvinBallOptimization.cpp.o'
    compile_command = ['/usr/bin/c++', *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'], *flags['CXX_FLAGS'],
                       '-o', str(obj), '-c', str(a.source/optimization)]
    subprocess.run(compile_command, check=True)
    link = shlex.split((cwd/'CMakeFiles/KelvinBall.dir/link.txt').read_text())
    if sum(token.endswith('KelvinBallOptimization.cpp.o') for token in link) != 1:
        raise RuntimeError('Unexpected original link command')
    link = [str(obj) if token.endswith('KelvinBallOptimization.cpp.o') else token for token in link]
    link[link.index('-o')+1] = str(a.build/'bin/KelvinBall')
    inputs = {}
    for token in link[1:]:
        path = (cwd/token).resolve()
        if path.is_file() and token.endswith(('.a','.so','.o')) and path != obj:
            inputs[str(path)] = sha(path)
    subprocess.run(link, cwd=cwd, check=True)
    if any(sha(Path(path)) != digest for path,digest in inputs.items()):
        raise RuntimeError('Immutable framework inputs changed during linking')
    shutil.copy2(oldbuild/'bin/KelvinCheckpointTest', a.build/'bin/KelvinCheckpointTest')
    write(a.build/'reused-framework-provenance.json', dict(
        source=str(a.source), commit=git(a.source,'rev-parse','HEAD'), framework_campaign=str(a.campaign),
        framework_source=str(oldsource), framework_commit=old['commit'], framework_build=str(oldbuild),
        compiler_command=compile_command, link_command=link, link_cwd=str(cwd),
        immutable_input_hashes=inputs, compiler_exit_code=0, linker_exit_code=0,
        binary_sha256=sha(a.build/'bin/KelvinBall'),
        checkpoint_test_sha256=sha(a.build/'bin/KelvinCheckpointTest')))


if __name__ == '__main__':
    main()
