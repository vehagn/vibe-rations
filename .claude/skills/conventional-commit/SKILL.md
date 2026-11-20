---
name: conventional-commit
description: Commit work in vibe-rations the way the maintainer wants - conventional commits, unsigned, one logical change each, on a branch, with `just ci` passing first. Use whenever you commit, and commit as you go rather than at the end.
---

# Committing

The maintainer wants commits made as work progresses, each a single logical change, in
the conventional-commits style, and **unsigned**: their global git config signs commits,
so pass `--no-gpg-sign` every time. Don't change their git config.

## Steps

1. Work on a branch, not `main` (`git switch -c <type>/<topic>`, e.g. `feat/adjoint`).
2. `just in ci` must pass: formatting, clippy with `-D warnings`, and every test.
3. Stage only the files that belong to this change, by path. Never stage `out/`, scratch
   tests or anything git-ignored.
4. Write the message to a file and commit from it. Heredocs inside `&&` chains have eaten
   subject lines before, so don't use them for messages:

   ```bash
   printf '%s\n' 'feat: back-propagate adjoint sources' '' \
     'Inject residuals at the receivers as forces, in reversed time.' '' \
     'Co-Authored-By: <the agent attribution line, if your instructions require one>' \
     > /tmp/msg.txt
   git commit --no-gpg-sign -F /tmp/msg.txt
   git log --format='%h %G? %s' -3     # %G? must be N: unsigned
   ```

5. Don't push unless asked.

## The message

```
<type>(<optional scope>): <imperative summary, lower case, no full stop, ≤ 72 chars>

<why, and what a reviewer needs to know, wrapped at 72>
```

| Type | For |
|---|---|
| `feat` | new behaviour: a solver feature, a configuration option, a CLI command |
| `fix` | a bug fix, including a wrong physical result |
| `perf` | faster or leaner, same results |
| `refactor` | restructuring, same behaviour |
| `test` | tests only |
| `docs` | README, AGENTS.md, skills, comments |
| `build` | Cargo, the dev container, the Justfile |
| `chore` | anything else |

Scopes, when they help: `solver`, `cpml`, `free-surface`, `source`, `receiver`, `config`,
`cli`, `fwi`. Breaking changes (configuration keys, output formats) get a `!` after the
type and a `BREAKING CHANGE:` footer.

If a physics change moves a baseline, say so in the body with the before and after
numbers.
