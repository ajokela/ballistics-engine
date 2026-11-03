//! The commit guard in `.githooks/` refuses what must not reach this public repository.
//!
//! Each test builds a throwaway repository under the system temporary directory, copies the
//! guard (`.githooks/` and `.sensitive-patterns`) into it, points `core.hooksPath` at the copy
//! and commits through the real `git`, so what is tested is the guard as git runs it.
//!
//! The local pattern list read by rule 3 of `pre-commit` and by `commit-msg` is never part of a
//! repository: it lives in the clone's git directory (`info/private-patterns`) or wherever
//! `hooks.privatePatterns` points. These tests write DUMMY patterns only, never a real name.
//!
//! Every test skips, saying so, when `git` or `bash` is not installed.
#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

/// The dummy name the patterns below stand for.
const TOKEN: &str = "PRIVATE_DUMMY_TOKEN";

/// The pattern in [`PATTERNS`] that matches [`TOKEN`], as the refusal quotes it.
const PATTERN: &str = "private[_-]?dummy[_-]?token";

/// A pattern list in the documented format. The comment, the empty line and the line of spaces
/// are ignored (read as patterns, they would match the ordinary text in [`ORDINARY`]), and so is
/// the whitespace around a pattern (kept, it would stop the patterns matching the test lines).
const PATTERNS: &str = "# ordinary\n\n   \n  private[_-]?dummy[_-]?token\nSECOND_DUMMY_NAME \t\n";

/// Text no pattern in [`PATTERNS`] matches.
const ORDINARY: &str = "# ordinary heading\nwords   spaced   out\n";

/// Environment a test run may inherit that would point `git` at another repository, add
/// configuration, or turn the guard off.
const INHERITED: [&str; 13] = [
    "GIT_DIR",
    "GIT_WORK_TREE",
    "GIT_INDEX_FILE",
    "GIT_OBJECT_DIRECTORY",
    "GIT_ALTERNATE_OBJECT_DIRECTORIES",
    "GIT_COMMON_DIR",
    "GIT_NAMESPACE",
    "GIT_PREFIX",
    "GIT_CONFIG",
    "GIT_CONFIG_PARAMETERS",
    "GIT_CONFIG_COUNT",
    "GIT_CEILING_DIRECTORIES",
    "SENSITIVE_OK",
];

static NEXT: AtomicUsize = AtomicUsize::new(0);

/// A directory under the system temporary directory, removed when dropped, so also when a test
/// panics.
struct Sandbox {
    root: PathBuf,
}

impl Sandbox {
    /// A new sandbox, or `None` (having said why) when the guard cannot be run here.
    fn new() -> Option<Sandbox> {
        let runs = |program: &str| {
            Command::new(program)
                .arg("--version")
                .output()
                .is_ok_and(|out| out.status.success())
        };
        if !(runs("git") && runs("bash")) {
            eprintln!("skipping the commit guard test: `git` or `bash` is not available");
            return None;
        }
        if !source().join(".githooks/pre-commit").is_file() {
            eprintln!("skipping the commit guard test: this source tree has no .githooks/");
            return None;
        }
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_or(0, |elapsed| elapsed.as_nanos());
        let root = std::env::temp_dir().join(format!(
            "commit-guard-{}-{}-{nanos}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let sandbox = Sandbox { root };
        let home = sandbox.home();
        fs::create_dir_all(&home).unwrap();
        // An empty global configuration, so the owner's own git settings do not reach the test.
        fs::write(home.join(".gitconfig"), "").unwrap();
        Some(sandbox)
    }

    fn home(&self) -> PathBuf {
        self.root.join("home")
    }

    /// A repository with the guard copied in, active and committed, and no pattern file yet.
    fn guarded_repo(&self) -> Checkout<'_> {
        let dir = self.root.join("repo");
        fs::create_dir_all(&dir).unwrap();
        let repo = Checkout { sandbox: self, dir };
        repo.ok(&["init", "-q"]);
        for (key, value) in [
            ("user.name", "Guard Test"),
            ("user.email", "guard-test@example.invalid"),
            ("commit.gpgsign", "false"),
            ("core.autocrlf", "false"),
            ("core.hooksPath", ".githooks"),
        ] {
            repo.ok(&["config", key, value]);
        }
        let hooks = repo.dir.join(".githooks");
        fs::create_dir_all(&hooks).unwrap();
        for entry in fs::read_dir(source().join(".githooks")).unwrap() {
            let entry = entry.unwrap();
            if entry.file_type().unwrap().is_file() {
                let copy = hooks.join(entry.file_name());
                fs::copy(entry.path(), &copy).unwrap();
                fs::set_permissions(&copy, fs::Permissions::from_mode(0o755)).unwrap();
            }
        }
        fs::copy(
            source().join(".sensitive-patterns"),
            repo.dir.join(".sensitive-patterns"),
        )
        .unwrap();
        repo.ok(&["add", "--", ".githooks", ".sensitive-patterns"]);
        repo.commit("Add the commit guard").expect_landed();
        repo
    }
}

impl Drop for Sandbox {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn source() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

/// One working tree of a sandboxed repository.
struct Checkout<'a> {
    sandbox: &'a Sandbox,
    dir: PathBuf,
}

/// The outcome of one `git commit`.
struct Attempt {
    succeeded: bool,
    landed: bool,
    output: String,
}

impl Attempt {
    /// The commit was refused by `hook`, whose messages say each of `says`.
    #[track_caller]
    fn expect_refused(&self, hook: &str, says: &[&str]) {
        assert!(
            !self.succeeded && !self.landed,
            "expected the guard to refuse this commit, but it went through:\n{}",
            self.output
        );
        let refusal = format!("{hook} GUARD: refusing this commit");
        for needle in says.iter().chain(&[refusal.as_str()]) {
            assert!(
                self.output.contains(needle),
                "the refusal does not say `{needle}`:\n{}",
                self.output
            );
        }
    }

    #[track_caller]
    fn expect_landed(&self) {
        assert!(
            self.succeeded && self.landed,
            "expected this commit to land:\n{}",
            self.output
        );
    }
}

fn text(out: &Output) -> String {
    format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    )
}

impl Checkout<'_> {
    /// `program` run in this checkout, with the sandbox's home and none of the inherited git
    /// state.
    fn command(&self, program: &str) -> Command {
        let mut command = Command::new(program);
        command.current_dir(&self.dir);
        for name in INHERITED {
            command.env_remove(name);
        }
        let home = self.sandbox.home();
        command
            .env("HOME", &home)
            .env("XDG_CONFIG_HOME", home.join(".config"))
            .env("GIT_CONFIG_GLOBAL", home.join(".gitconfig"))
            .env("GIT_CONFIG_NOSYSTEM", "1")
            .env("GIT_EDITOR", ":")
            .env("GIT_TERMINAL_PROMPT", "0");
        command
    }

    fn git(&self, args: &[&str], env: &[(&str, &str)]) -> Output {
        let mut command = self.command("git");
        command.args(args);
        for (name, value) in env {
            command.env(name, value);
        }
        command
            .output()
            .unwrap_or_else(|error| panic!("git {args:?}: {error}"))
    }

    #[track_caller]
    fn ok(&self, args: &[&str]) -> String {
        let out = self.git(args, &[]);
        assert!(out.status.success(), "git {args:?} failed:\n{}", text(&out));
        String::from_utf8_lossy(&out.stdout).into_owned()
    }

    fn write(&self, path: &str, contents: impl AsRef<[u8]>) {
        let path = self.dir.join(path);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, contents).unwrap();
    }

    fn stage(&self, path: &str, contents: impl AsRef<[u8]>) {
        self.write(path, contents);
        self.ok(&["add", "--", path]);
    }

    fn head(&self) -> String {
        let out = self.git(&["rev-parse", "-q", "--verify", "HEAD"], &[]);
        String::from_utf8_lossy(&out.stdout).trim().to_string()
    }

    fn attempt(&self, args: &[&str], env: &[(&str, &str)]) -> Attempt {
        let before = self.head();
        let out = self.git(args, env);
        Attempt {
            succeeded: out.status.success(),
            landed: self.head() != before,
            output: text(&out),
        }
    }

    fn commit(&self, message: &str) -> Attempt {
        self.attempt(&["commit", "-q", "-m", message], &[])
    }

    /// Commit `contents` as `path` past the guard, as a reviewed exception would.
    fn commit_past_the_guard(&self, path: &str, contents: &str) {
        self.stage(path, contents);
        self.attempt(
            &["commit", "-q", "-m", "A reviewed exception"],
            &[("SENSITIVE_OK", "1")],
        )
        .expect_landed();
    }

    /// Throw away whatever a refused commit left staged.
    fn discard(&self) {
        self.ok(&["reset", "-q", "--hard"]);
        self.ok(&["clean", "-fdq"]);
    }

    /// Write the pattern file at its default place: `info/private-patterns` in the clone's
    /// common git directory, which every worktree of the clone shares.
    fn private_patterns(&self, patterns: &str) {
        let common = self.ok(&["rev-parse", "--git-common-dir"]);
        let info = self.dir.join(common.trim()).join("info");
        fs::create_dir_all(&info).unwrap();
        fs::write(info.join("private-patterns"), patterns).unwrap();
    }

    fn last_message(&self) -> String {
        self.ok(&["log", "-1", "--format=%B"]).trim().to_string()
    }
}

/// `line 1\n` through `line N\n`.
fn numbered(lines: std::ops::RangeInclusive<u32>) -> String {
    lines.map(|n| format!("line {n}\n")).collect()
}

/// A UTF-8 locale this machine has, or `None` (having said so) when it has none. A committer
/// usually runs in one, and in one the BSD `grep` of macOS skips every line that is not valid
/// UTF-8, where GNU `grep` reads it.
fn utf8_locale() -> Option<String> {
    let listed = Command::new("locale")
        .arg("-a")
        .output()
        .map(|out| String::from_utf8_lossy(&out.stdout).into_owned())
        .unwrap_or_default();
    let found = ["en_US.UTF-8", "C.UTF-8", "en_US.utf8", "C.utf8"]
        .into_iter()
        .find(|name| listed.lines().any(|line| line.trim() == *name));
    if found.is_none() {
        eprintln!("no UTF-8 locale here; this test runs in the inherited locale instead");
    }
    found.map(str::to_string)
}

#[test]
fn an_added_line_matching_a_private_pattern_is_refused_naming_its_file_line_and_pattern() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    repo.stage("docs/notes.txt", numbered(1..=6));
    repo.stage("docs/moved.txt", numbered(1..=10));
    repo.commit("Add two ordinary files").expect_landed();

    // A modified file, where the matching line is the second of two added lines and line 7 of
    // the staged file; a new file; and a renamed file with a line added at its end (git must see
    // the rename for this to test one).
    repo.write(
        "docs/notes.txt",
        format!(
            "{}an ordinary new line\nthe Private-Dummy_token came in here\nline 6\n",
            numbered(1..=5)
        ),
    );
    repo.write("docs/new.txt", format!("{TOKEN}\n"));
    repo.ok(&["mv", "docs/moved.txt", "docs/renamed.txt"]);
    repo.write(
        "docs/renamed.txt",
        format!("{}and SECOND_DUMMY_NAME too\n", numbered(1..=10)),
    );
    repo.ok(&["add", "--", "docs"]);
    let status = repo.ok(&["diff", "--cached", "--name-status", "-M"]);
    assert!(
        status.contains("docs/moved.txt\tdocs/renamed.txt"),
        "git did not record a rename:\n{status}"
    );
    repo.commit("Add more to the files").expect_refused(
        "pre-commit",
        &[
            "docs/notes.txt:7 ",
            "docs/new.txt:1 ",
            "docs/renamed.txt:11 ",
            PATTERN,
            "SECOND_DUMMY_NAME",
        ],
    );
}

#[test]
fn a_binary_file_is_skipped() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    // The NUL bytes make git call the file binary; the matching line holds none of them, so
    // the file is passed only because it is skipped.
    let mut bytes = format!("{TOKEN}\n").into_bytes();
    bytes.extend_from_slice(&[0, 1, 2, 0]);
    repo.stage("data/blob.bin", bytes);
    repo.commit("Add a binary blob").expect_landed();
}

#[test]
fn a_commit_message_matching_a_private_pattern_is_refused() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    repo.stage("docs/ordinary.txt", ORDINARY);
    repo.commit("Add ordinary words\n\nThis names a Private-Dummy-Token.")
        .expect_refused("commit-msg", &["commit message line 3 ", PATTERN]);
    // `-m` keeps a `#` line (git strips comment lines from an edited message only), so it is read.
    repo.commit("Add ordinary words\n\n# SECOND_DUMMY_NAME")
        .expect_refused(
            "commit-msg",
            &["commit message line 3 ", "SECOND_DUMMY_NAME"],
        );
    repo.commit("Add ordinary words").expect_landed();
}

#[test]
fn comment_lines_that_git_strips_from_the_message_are_not_read() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    repo.commit_past_the_guard("docs/old.txt", &format!("keep this\n{TOKEN} goes\n"));
    repo.stage("docs/old.txt", "keep this\n");

    // An edited, verbose commit: git strips the comment line the editor adds and the diff below
    // the scissors line, which here removes a line matching a pattern.
    let editor = sandbox.root.join("editor.sh");
    fs::write(
        &editor,
        format!(
            "#!/bin/sh\n{{ printf 'Remove a line\\n# a comment naming {TOKEN}\\n'; cat \"$1\"; }} \
             > \"$1.edited\" && mv \"$1.edited\" \"$1\"\n"
        ),
    )
    .unwrap();
    fs::set_permissions(&editor, fs::Permissions::from_mode(0o755)).unwrap();
    repo.attempt(
        &["commit", "-q", "-v"],
        &[("GIT_EDITOR", editor.to_str().unwrap())],
    )
    .expect_landed();
    assert_eq!(repo.last_message(), "Remove a line");

    // With commit.cleanup=strip git strips comment lines from a `-m` message too.
    repo.ok(&["config", "commit.cleanup", "strip"]);
    repo.stage("docs/more.txt", ORDINARY);
    repo.commit("Add more\n# SECOND_DUMMY_NAME in a stripped comment")
        .expect_landed();
    assert_eq!(repo.last_message(), "Add more");
}

#[test]
fn an_unrelated_line_and_message_commit() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    repo.stage("docs/ordinary.txt", ORDINARY);
    repo.commit("Add ordinary words\n\n# ordinary heading, words   spaced   out")
        .expect_landed();
}

#[test]
fn a_line_already_in_the_file_is_not_read_again_until_it_changes() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    repo.commit_past_the_guard("docs/old.txt", &format!("{TOKEN} was here\nsecond line\n"));

    // Kept while an ordinary line is added below it.
    repo.stage(
        "docs/old.txt",
        format!("{TOKEN} was here\nsecond line\nthird line\n"),
    );
    repo.commit("Add a third line").expect_landed();

    // Carried along unchanged by a rename.
    repo.ok(&["mv", "docs/old.txt", "docs/kept.txt"]);
    repo.commit("Rename the file").expect_landed();

    // An edited line is an added line.
    repo.stage(
        "docs/kept.txt",
        format!("{TOKEN} was here, edited\nsecond line\nthird line\n"),
    );
    repo.commit("Edit the first line")
        .expect_refused("pre-commit", &["docs/kept.txt:1 ", PATTERN]);
    repo.discard();

    // A removed line is not added text.
    repo.stage("docs/kept.txt", "second line\nthird line\n");
    repo.commit("Remove the first line").expect_landed();
}

#[test]
fn without_a_pattern_file_rule_three_is_off() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.stage("docs/notes.txt", format!("{TOKEN}\n"));
    repo.commit(&format!("Name {TOKEN}")).expect_landed();

    // A configured path that does not exist turns the rule off as well, and says so.
    let missing = sandbox.root.join("no-such-patterns");
    let missing = missing.to_str().unwrap();
    repo.ok(&["config", "hooks.privatePatterns", missing]);
    repo.stage("docs/more.txt", format!("{TOKEN} again\n"));
    let attempt = repo.commit(&format!("Name {TOKEN} again"));
    attempt.expect_landed();
    assert!(
        attempt
            .output
            .contains(&format!("'{missing}', which does not exist")),
        "the hook did not say the configured pattern file is missing:\n{}",
        attempt.output
    );
}

#[test]
fn hooks_private_patterns_overrides_the_default_path() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns("DEFAULT_ONLY_DUMMY\n");
    let elsewhere = sandbox.root.join("elsewhere-patterns");
    fs::write(&elsewhere, "OVERRIDE_ONLY_DUMMY\n").unwrap();
    repo.ok(&[
        "config",
        "hooks.privatePatterns",
        elsewhere.to_str().unwrap(),
    ]);

    repo.stage("docs/a.txt", "DEFAULT_ONLY_DUMMY\n");
    repo.commit("Name DEFAULT_ONLY_DUMMY").expect_landed();
    repo.stage("docs/b.txt", "OVERRIDE_ONLY_DUMMY\n");
    repo.commit("Add b")
        .expect_refused("pre-commit", &["docs/b.txt:1 ", "OVERRIDE_ONLY_DUMMY"]);
    repo.discard();
    repo.stage("docs/c.txt", ORDINARY);
    repo.commit("Name OVERRIDE_ONLY_DUMMY").expect_refused(
        "commit-msg",
        &["commit message line 1 ", "OVERRIDE_ONLY_DUMMY"],
    );

    // A `~/` path is the user's home directory, as for any git path setting.
    fs::write(sandbox.home().join("patterns-at-home"), "HOME_DUMMY\n").unwrap();
    repo.ok(&["config", "hooks.privatePatterns", "~/patterns-at-home"]);
    repo.commit("Name HOME_DUMMY")
        .expect_refused("commit-msg", &["commit message line 1 ", "HOME_DUMMY"]);
    repo.commit("Name OVERRIDE_ONLY_DUMMY").expect_landed();
}

#[test]
fn sensitive_ok_lets_a_reviewed_exception_through_both_hooks() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    repo.stage("docs/notes.txt", format!("{TOKEN}\n"));
    let message = format!("Name {TOKEN}");
    repo.commit(&message)
        .expect_refused("pre-commit", &[PATTERN]);
    repo.attempt(&["commit", "-q", "-m", &message], &[("SENSITIVE_OK", "1")])
        .expect_landed();
}

#[test]
fn the_earlier_rules_still_refuse() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);

    // A staged path matching .sensitive-patterns.
    repo.stage("notes/hornady_loads.txt", ORDINARY);
    repo.commit("Add loads").expect_refused(
        "pre-commit",
        &["notes/hornady_loads.txt", "matches sensitive pattern"],
    );
    repo.discard();

    // A new file naming the data's owner (spelled apart here, so this file does not).
    let party = ["Applied", "Ballistics"].join(" ");
    repo.stage("docs/source.md", format!("Data from {party}.\n"));
    repo.commit("Add a source")
        .expect_refused("pre-commit", &["docs/source.md", &format!("names {party}")]);
    repo.discard();

    // A new file about the patent-pending solver (spelled apart for the same reason).
    let solver = ["BC", "5D"].concat();
    repo.stage("docs/table.txt", format!("{solver} table\n"));
    repo.commit("Add a table").expect_refused(
        "pre-commit",
        &["docs/table.txt", &format!("contains {solver} material")],
    );
}

#[test]
fn the_default_pattern_file_covers_every_worktree_of_the_clone() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    let other = sandbox.root.join("other-worktree");
    repo.ok(&[
        "worktree",
        "add",
        "-q",
        "-b",
        "side",
        other.to_str().unwrap(),
    ]);
    let other = Checkout {
        sandbox: &sandbox,
        dir: other,
    };
    other.stage("docs/notes.txt", format!("{TOKEN}\n"));
    other
        .commit("Add notes")
        .expect_refused("pre-commit", &["docs/notes.txt:1 ", PATTERN]);
}

#[test]
fn a_copy_of_the_pattern_file_cannot_be_staged() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    // A regex does not match its own text, so only the copy check can catch this list.
    let patterns = "# ordinary\nprivate[_-]?dummy[_-]?token\n";
    repo.private_patterns(patterns);
    repo.stage("docs/patterns.txt", patterns);
    repo.commit("Keep a copy of the list").expect_refused(
        "pre-commit",
        &["docs/patterns.txt", "copy of the private pattern file"],
    );
}

#[test]
fn a_pattern_that_is_not_a_valid_regex_refuses_rather_than_passing() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns("(unbalanced\n");
    repo.stage("docs/ordinary.txt", ORDINARY);
    repo.commit("Add ordinary words").expect_refused(
        "pre-commit",
        &["(unbalanced", "not a valid extended regular expression"],
    );

    // commit-msg on its own, as git would run it for a `-m` message.
    let message = sandbox.root.join("message.txt");
    fs::write(&message, "Add ordinary words\n").unwrap();
    let out = repo
        .command("bash")
        .arg(".githooks/commit-msg")
        .arg(&message)
        .output()
        .unwrap();
    let output = text(&out);
    assert!(
        !out.status.success() && output.contains("not a valid extended regular expression"),
        "commit-msg did not refuse an invalid pattern:\n{output}"
    );
}

#[test]
fn a_matching_line_that_is_not_valid_utf8_is_still_refused_in_a_utf8_locale() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    let locale = utf8_locale();
    let env: Vec<(&str, &str)> = locale
        .as_deref()
        .map(|name| vec![("LC_ALL", name), ("LANG", name)])
        .unwrap_or_default();

    // Latin-1 text: 0xE9 is "é" there and no valid UTF-8 sequence.
    repo.stage(
        "docs/latin1.txt",
        [
            &b"caf\xe9 au lait\n"[..],
            b"notes caf\xe9 ",
            TOKEN.as_bytes(),
            b"\n",
        ]
        .concat(),
    );
    repo.attempt(&["commit", "-q", "-m", "Add notes"], &env)
        .expect_refused("pre-commit", &["docs/latin1.txt:2 ", PATTERN]);
    repo.discard();

    // The same text in a message given with -F.
    repo.stage("docs/ordinary.txt", ORDINARY);
    let message = sandbox.root.join("latin1-message.txt");
    fs::write(
        &message,
        [&b"Add notes\n\nnotes caf\xe9 "[..], TOKEN.as_bytes(), b"\n"].concat(),
    )
    .unwrap();
    repo.attempt(&["commit", "-q", "-F", message.to_str().unwrap()], &env)
        .expect_refused("commit-msg", &["commit message line 3 ", PATTERN]);
}

#[test]
fn an_exported_cdpath_does_not_move_the_lookup_of_the_default_pattern_file() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    // In the main checkout git names the common git directory by the relative `.git`, which `cd`
    // would look up through CDPATH: here first in a decoy directory holding a `.git` of its own,
    // whose list names something else.
    let decoy = sandbox.root.join("decoy");
    fs::create_dir_all(decoy.join(".git/info")).unwrap();
    fs::write(
        decoy.join(".git/info/private-patterns"),
        "DECOY_ONLY_DUMMY\n",
    )
    .unwrap();
    let cdpath = format!("{}:.", decoy.to_str().unwrap());
    let env = [("CDPATH", cdpath.as_str())];

    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    let common = repo.ok(&["rev-parse", "--git-common-dir"]);
    assert_eq!(
        common.trim(),
        ".git",
        "this test needs a relative git directory"
    );

    repo.stage("docs/notes.txt", format!("{TOKEN}\n"));
    repo.attempt(&["commit", "-q", "-m", "Add notes"], &env)
        .expect_refused("pre-commit", &["docs/notes.txt:1 ", PATTERN]);
    repo.discard();

    repo.stage("docs/ordinary.txt", ORDINARY);
    let message = format!("Name {TOKEN}");
    repo.attempt(&["commit", "-q", "-m", &message], &env)
        .expect_refused("commit-msg", &["commit message line 1 ", PATTERN]);
    // Only the clone's own list is read.
    repo.attempt(&["commit", "-q", "-m", "Name DECOY_ONLY_DUMMY"], &env)
        .expect_landed();
}

#[test]
fn grep_options_from_the_environment_do_not_change_what_matches() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    // The BSD grep of macOS reads options from GREP_OPTIONS; -x would let only a whole-line
    // match through.
    let env = [("GREP_OPTIONS", "-x")];
    repo.stage("docs/notes.txt", format!("a line naming {TOKEN}\n"));
    repo.attempt(&["commit", "-q", "-m", "Add notes"], &env)
        .expect_refused("pre-commit", &["docs/notes.txt:1 ", PATTERN]);
    repo.discard();

    repo.stage("docs/ordinary.txt", ORDINARY);
    let message = format!("Add words\n\nThey name {TOKEN}.");
    repo.attempt(&["commit", "-q", "-m", &message], &env)
        .expect_refused("commit-msg", &["commit message line 3 ", PATTERN]);
}

#[test]
fn the_line_reported_is_right_when_the_diff_carries_context_lines() {
    let Some(sandbox) = Sandbox::new() else {
        return;
    };
    let repo = sandbox.guarded_repo();
    repo.private_patterns(PATTERNS);
    // Thirty numbered lines, line 18 blank.
    let mut lines: Vec<String> = numbered(1..=30).lines().map(str::to_string).collect();
    lines[17] = String::new();
    repo.stage("docs/counted.txt", lines.join("\n") + "\n");
    repo.commit("Add a counted file").expect_landed();

    // Lines 10 and 20 change; only line 20 matches. diff.interHunkContext joins the two hunks
    // across the context between them, GIT_DIFF_OPTS overrides -U0 with context of its own, and
    // diff.suppressBlankEmpty prints the blank context line 18 as an empty line.
    lines[9] = "an ordinary tenth line".to_string();
    lines[19] = format!("the {TOKEN} on line twenty");
    repo.stage("docs/counted.txt", lines.join("\n") + "\n");
    let says = ["docs/counted.txt:20 ", PATTERN];
    let more_context = [("GIT_DIFF_OPTS", "--unified=3")];

    repo.ok(&["config", "diff.interHunkContext", "10"]);
    repo.commit("Edit two lines")
        .expect_refused("pre-commit", &says);
    repo.ok(&["config", "--unset", "diff.interHunkContext"]);
    repo.attempt(&["commit", "-q", "-m", "Edit two lines"], &more_context)
        .expect_refused("pre-commit", &says);
    repo.ok(&["config", "diff.suppressBlankEmpty", "true"]);
    repo.attempt(&["commit", "-q", "-m", "Edit two lines"], &more_context)
        .expect_refused("pre-commit", &says);
}

#[test]
fn both_hooks_carry_the_same_shared_block() {
    // The pattern file's lookup, its reading and the matching are one block of text in both
    // hooks, so a fix to one cannot leave the other behind.
    let hooks = source().join(".githooks");
    if !hooks.join("pre-commit").is_file() {
        eprintln!("skipping the commit guard test: this source tree has no .githooks/");
        return;
    }
    let block = |hook: &str| {
        let text = fs::read_to_string(hooks.join(hook)).unwrap();
        let start = text
            .find("# --- BEGIN shared block")
            .unwrap_or_else(|| panic!("{hook} has no shared block"));
        let end = text[start..]
            .find("# --- END shared block")
            .unwrap_or_else(|| panic!("{hook}'s shared block has no end"));
        text[start..start + end].to_string()
    };
    assert_eq!(
        block("pre-commit"),
        block("commit-msg"),
        "the shared block differs between .githooks/pre-commit and .githooks/commit-msg"
    );
}
