# Known Issues and Fixes

## Garden Matplotlib Compatibility Issues

**Problems:** 
Garden matplotlib is incompatible with newer matplotlib versions. Missing event handler methods:
- `resize_event`
- `motion_notify_event`
- `button_press_event`
- `button_release_event`
- `scroll_event`

**Solution:** These methods are automatically disabled during setup. The visualization still works but without interactive features (zoom, pan, etc).

**Manual patch** (if needed after reinstalling garden matplotlib):
```bash
cd /path/to/miedema
./setup.sh
```

Or run the patch manually:
```python
import os

garden_path = os.path.expanduser("~/.kivy/garden/garden.matplotlib/backend_kivy.py")

with open(garden_path, 'r') as f:
    content = f.read()

lines = content.split('\n')
fixed_lines = []
skip_next = False

for i, line in enumerate(lines):
    if skip_next:
        skip_next = False
        continue
    
    if 'self.resize_event()' in line and 'def resize_event' not in line:
        indent = len(line) - len(line.lstrip())
        fixed_lines.append(' ' * indent + 'pass  # disabled')
    elif 'self.motion_notify_event(' in line and 'def motion_notify_event' not in line:
        indent = len(line) - len(line.lstrip())
        fixed_lines.append(' ' * indent + 'pass  # disabled')
    elif 'self.button_press_event(' in line and 'def button_press_event' not in line:
        indent = len(line) - len(line.lstrip())
        fixed_lines.append(' ' * indent + 'pass  # disabled')
        if i + 1 < len(lines) and 'dblclick' in lines[i + 1]:
            skip_next = True
    elif 'self.button_release_event(' in line and 'def button_release_event' not in line:
        indent = len(line) - len(line.lstrip())
        fixed_lines.append(' ' * indent + 'pass  # disabled')
        if i + 1 < len(lines) and 'guiEvent' in lines[i + 1]:
            skip_next = True
    elif 'self.scroll_event(' in line and 'def scroll_event' not in line:
        indent = len(line) - len(line.lstrip())
        fixed_lines.append(' ' * indent + 'pass  # disabled')
    else:
        fixed_lines.append(line)

with open(garden_path, 'w') as f:
    f.write('\n'.join(fixed_lines))
```

**Alternative:** Consider using a different visualization backend in future versions.

## Cutoff Calculation

The automatic cutoff calculation uses `min(cell_params) * 0.8` instead of volume-based calculation for better neighbor detection.

Use `-1` as cutoff value to enable automatic calculation, or specify a custom value (e.g., `2.0`, `3.5`).

