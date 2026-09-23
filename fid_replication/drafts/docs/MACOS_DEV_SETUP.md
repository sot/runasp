# macOS development setup

The `ds` repo contains files whose names differ only in case
(e.g. `Aca_Id_Image.cc` and `aca_id_image.cc` in the same directory).
APFS is case-insensitive by default: if you check this repo out onto
a stock macOS volume, the second filename in each pair silently
overwrites the first, and one file is lost. Sometimes the lost file
is a stale duplicate; sometimes it is the live implementation.

This document is the breadcrumb trail for setting up a Mac so you
do not lose source files when working with this repo.

## 1. Detect collisions

`scripts/find_case_collisions.py` reads a `find` listing of the
upstream source and reports every group of names that collide on a
case-insensitive volume.

To regenerate the inventory:

```sh
# from a case-sensitive checkout of the repo
( cd /Volumes/chandra-ds-source/ds && find * ) > ls_all.txt
python3 scripts/find_case_collisions.py ls_all.txt > scripts/case_collisions.txt
```

`scripts/case_collisions.txt` (current snapshot, 2026-04-30) records
65 groups: 24 source-side and 41 in `build/`. The full list of
source-side collisions is in that file.

## 2. The fix: a case-sensitive sparse bundle

A separate APFS *volume* in the system disk container would be
cleaner, but on Apple Silicon you typically hit
`Error: -69493: You can't add any more APFS Volumes to its APFS Container`.
A sparse bundle is the workaround — a self-contained file that
behaves like a volume, grows on demand, and bypasses the container
limit.

```sh
hdiutil create -size 50g -fs "Case-sensitive APFS" \
    -type SPARSEBUNDLE -volname chandra-ds-source ~/chandra-ds-source
```

This creates `~/chandra-ds-source.sparsebundle`. Mount it manually with:

```sh
hdiutil attach ~/chandra-ds-source.sparsebundle    # mounts at /Volumes/chandra-ds-source
hdiutil detach /Volumes/chandra-ds-source          # unmount
```

Clone the repo into `/Volumes/chandra-ds-source/ds` and work there.
All file operations on that volume are case-sensitive, so the
upstream sources coexist correctly.

The `~/SAO/data_systems/ds` checkout on the case-insensitive home
volume is the *audit* copy — useful for confirming what's missing,
not for production work.

## 3. Auto-mount at login

`~/Library/LaunchAgents/com.javierg.chandra-ds-source.plist` runs
`hdiutil attach` at login so the volume is always available:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>Label</key>
    <string>com.javierg.chandra-ds-source</string>
    <key>ProgramArguments</key>
    <array>
        <string>/usr/bin/hdiutil</string>
        <string>attach</string>
        <string>-nobrowse</string>
        <string>/Users/javierg/chandra-ds-source.sparsebundle</string>
    </array>
    <key>RunAtLoad</key>
    <true/>
    <key>StandardOutPath</key>
    <string>/Users/javierg/Library/Logs/chandra-ds-source-mount.log</string>
    <key>StandardErrorPath</key>
    <string>/Users/javierg/Library/Logs/chandra-ds-source-mount.log</string>
</dict>
</plist>
```

To activate without logging out:

```sh
launchctl load ~/Library/LaunchAgents/com.javierg.chandra-ds-source.plist
```

To remove the auto-mount:

```sh
launchctl unload ~/Library/LaunchAgents/com.javierg.chandra-ds-source.plist
rm ~/Library/LaunchAgents/com.javierg.chandra-ds-source.plist
```

`-nobrowse` keeps the volume off the Finder sidebar. Drop that
argument if you want the volume visible.

## 4. Verify

```sh
mount | grep chandra-ds-source                          # should show /Volumes/chandra-ds-source
diskutil info /Volumes/chandra-ds-source | grep -i case # "Case-sensitive: Yes"
tail ~/Library/Logs/chandra-ds-source-mount.log         # check launchd output if mount failed
```

A quick functional check — these two paths should be distinct files
on the case-sensitive volume:

```sh
cd /Volumes/chandra-ds-source/ds
ls dstools/asp/aca_id_image/Aca_Id_Image.cc dstools/asp/aca_id_image/aca_id_image.cc
```

## 5. Recovery workflow when collisions strike

If you find a file missing on the case-insensitive copy:

1. Confirm it's a case-collision casualty: `grep <basename> scripts/case_collisions.txt`.
2. From the case-sensitive checkout at `/Volumes/chandra-ds-source/ds`,
   copy the surviving (mixed-case) file back into your audit working
   tree only if you genuinely need it there — usually you do not,
   since the case-sensitive volume is now your real workspace.
3. If you need to remove the local-volume casualties before
   re-cloning, `scripts/remove_case_collision_files.sh` lists every
   path. `rm -f` makes it safe even though APFS resolves both names
   to a single inode.

## 6. Why not just rename one file?

Renaming `Aca_Id_Image.cc` → `Aca_Id_Image_main.cc` would break the
collision but creates a permanent local fork: each renamed file
forces edits to its `CMakeLists.txt` `SRCS` line and to every
`#include` that names it, and the patch must be kept out of upstream
commits forever. Survivable for a one-off audit, painful as an
ongoing workspace. Use the sparse bundle.
