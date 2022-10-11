# Contributing

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.

## Environment setup

Nothing easier!

Fork and clone the repository:

```bash
git clone https://github.com/lgmgeo/AnnotSV
cd AnnotSV
```

We use [pdm](https://pdm.fming.dev) to manage the project and dependencies, install PDM if it isn't done yet, then:

```bash
pdm install
```

You now have the dependencies installed.

You can run the tests with `pdm run test [ARGS...]`.

## Test against multiple Python versions

This project uses [nox](https://nox.thea.codes/) as the test runner. See what sessions are list:

```bash
nox --list
```

And run the test suite on specified Python versions:

```bash
nox -s tests-3.8
```

## Pull requests guidelines

Link to any related issue in the Pull Request message.

During review, we recommend using fixups:

```bash
# SHA is the SHA of the commit you want to fix
git commit --fixup=SHA
```

Once all the changes are approved, you can squash your commits:

```bash
git rebase -i --autosquash master
```

And force-push:

```bash
git push -f
```

If this seems all too complicated, you can push or force-push each new commit,
and we will squash them ourselves if needed, before merging.
