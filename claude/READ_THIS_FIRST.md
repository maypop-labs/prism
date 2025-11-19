# READ THIS FIRST: General MCP Best Practices

**Audience:** Claude (you, reading this in the future)  
**Purpose:** Core lessons learned across multiple projects to avoid repeating mistakes  
**Last Updated:** October 28, 2025

---

## Critical Lessons Summary

1. **Token Conservation** - Don't exhaust context windows through careless bulk reads
2. **MCP Filesystem Efficiency** - Use the right tools consistently, don't mix systems
3. **Incremental Progress** - Work in small batches with user checkpoints

---

## Lesson 1: Token Conservation

### The Problem
In "Project file management overview" conversation (PRISM project), hit token limit on **second user prompt** because I read all 8 manager files simultaneously (10,000+ lines) and quoted large code blocks in my analysis.

### Core Principles

**1. Reconnaissance Before Reading**
- Always `list_directory` or `search_files` before reading files
- Check file sizes with `get_file_info` - if >500 lines, reconsider bulk reads
- Never read multiple large files simultaneously unless editing them that same turn

**2. Just Enough Information**
- Use `view` with line ranges for specific sections, not entire files
- Use `search_files` with patterns instead of reading entire files
- Summarize findings, don't quote code blocks
- Present counts/summaries, offer details on request

**3. Incremental Analysis**
- Analyze 1-3 files at a time, then checkpoint with user
- Don't analyze everything at once - user may redirect you
- Pattern: Analyze → Present summary → Wait for user → Continue

**4. Output Economy**
- Tables and lists beat prose
- Function names only, not bodies: `functionName(params)`
- Defer details: "Found 42 items - want specifics?" not a full list

**5. Red Flags (Stop and Reconsider)**
- About to read >3 files at once?
- About to quote >20 lines of code?
- Response draft is >500 words?
- Multiple tool calls in a row without pausing?

**6. Edit Mode Exception**
Only bulk-read multiple large files when you're about to edit them in that same turn.

### Quick Decision Tree
```
User asks question
    ↓
Can I answer from context?
    ├─ YES → Answer without tools
    └─ NO → Do I need file contents?
        ├─ NO → Use list/search/info only
        └─ YES → Do I need >1 file?
            ├─ NO → Read the one file
            └─ YES → How many?
                ├─ 2-3 small files → OK
                ├─ 2-3 large files → One at a time
                └─ >3 files → Batch in groups of 2-3
```

---

## Lesson 2: MCP Filesystem Efficiency

### The Problem
In SOLID project, took 5+ minutes and excessive tool calls to change ONE line of CSS because I was mixing tool systems (MCP Filesystem vs bash tools) and using wrong path formats.

### Core Principles

**1. Use MCP Filesystem Tools Consistently**
When you have MCP filesystem access, use MCP tools (either `Filesystem:` or `filesystem:` namespace):
- `list_directory` - List contents
- `read_file` - Read files
- `write_file` - Create new files
- `edit_file` - **Edit existing files** ← USE THIS
- `search_files` - Find files by pattern
- `search_file_contents` - Search within file contents (grep-like, available in lowercase namespace)
- `get_file_info` - Check file metadata
- Use whichever namespace has the most useful tools for your current task

**2. Don't Mix Tool Systems**
- ❌ Don't use bash `str_replace` when you have MCP filesystem access
- ❌ Don't use bash `grep` or `find` - use MCP `search_files` or `search_file_contents`
- ❌ Don't switch between bash and MCP tools for the same operation
- ✅ Pick one system (MCP Filesystem) and use it consistently
- ✅ Either `Filesystem:` or `filesystem:` namespace is fine - use the one with the tools you need

**3. Path Format Matters**
Always use forward slashes, even on Windows:
- ✅ `E:/lib/solid/styles/main.css`
- ❌ `E:\lib\solid\styles\main.css`
- ❌ `/mnt/project/styles/main.css` (Unix paths don't work in Windows MCP context)

**4. Efficient Edit Workflow**
To edit a file:
1. `read_file` (if you need to see current content)
2. `edit_file` with oldText/newText pairs
3. Done. That's it. Two tool calls maximum.

**5. Don't Overthink**
The direct approach is usually correct:
- Need to edit a file? Use `edit_file`
- Need to find files? Use `search_files`
- Need to search file contents? Use `search_file_contents` (lowercase namespace)
- Need to list a directory? Use `list_directory`

### Example: Editing a File (Right Way)

User: "Change the inline code color to red in markdown.css"

**Efficient (2 tool calls, 30 seconds):**
```
1. read_file → See current styling
2. edit_file → Change color
3. Done
```

**Inefficient (what not to do):**
```
1. Try bash grep
2. Try Unix paths that don't exist
3. Try str_replace with wrong format
4. Try Filesystem with backslashes
5. Finally succeed after 5 minutes
6. User is frustrated
```

---

## Lesson 3: General Efficiency Principles

### Pattern Recognition
Both token conservation and filesystem efficiency issues share root causes:
- **Mixing approaches** instead of using one consistent method
- **Overthinking** instead of using the direct solution
- **Bulk operations** instead of incremental work
- **Lack of planning** before executing tool calls

### Success Pattern
1. **Plan first** - What's the minimum information needed?
2. **Use one tool system** - Pick the right tool and stick with it
3. **Work incrementally** - Small batches with user checkpoints
4. **Be concise** - Tables/summaries beat verbose prose
5. **Watch red flags** - Multiple approaches? Bulk reads? Long responses?

### Universal Red Flags
- Switching between different tool approaches for same task
- Reading/processing more than you immediately need
- Response getting long without user checkpoint
- More than 3-4 tool calls in a row

---

## Quick Reference Card

**Before any major operation, ask yourself:**

1. □ Do I need file contents, or just metadata?
2. □ Am I using MCP Filesystem tools consistently?
3. □ Am I reading only what I need right now?
4. □ Have I checked file sizes before bulk reads?
5. □ Am I working incrementally with checkpoints?
6. □ Is my planned response concise and actionable?

**If you answered NO to any question, reconsider your approach.**

---

## Memory Anchors

**Token Conservation:** "Project file management overview" died on prompt 2 because I loaded 10,000+ lines before starting work.

**Filesystem Efficiency:** Took 5 minutes to change one CSS line because I mixed tool systems instead of using `edit_file` directly.

**Core Insight:** Both mistakes stem from the same problem - not using the right tool in the right way at the right scale.

---

## When in Doubt

- **Token Conservation:** Do less, ask "want me to continue?"
- **File Operations:** Use MCP `edit_file`, it just works
- **General Principle:** Simple, direct, incremental beats complex, bulk, all-at-once

**Remember:** Efficiency isn't about speed - it's about not wasting resources (tokens, time, user patience) through poor tool choices.
