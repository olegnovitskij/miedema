#!/bin/bash

set -e

echo "Creating virtual environment..."
uv venv

echo "Activating virtual environment..."
source .venv/bin/activate

echo "Installing dependencies..."
uv pip install -e .

echo "Installing kivy garden matplotlib..."
garden install matplotlib

echo "Applying garden.matplotlib patch..."
python << 'PATCH'
import os

garden_path = os.path.expanduser("~/.kivy/garden/garden.matplotlib/backend_kivy.py")

if os.path.exists(garden_path):
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
            fixed_lines.append(' ' * indent + 'pass  # self.resize_event()')
        elif 'self.motion_notify_event(' in line and 'def motion_notify_event' not in line:
            indent = len(line) - len(line.lstrip())
            fixed_lines.append(' ' * indent + 'pass  # self.motion_notify_event(...)')
        elif 'self.button_press_event(' in line and 'def button_press_event' not in line:
            indent = len(line) - len(line.lstrip())
            fixed_lines.append(' ' * indent + 'pass  # self.button_press_event(...)')
            if i + 1 < len(lines) and 'dblclick' in lines[i + 1]:
                skip_next = True
        elif 'self.button_release_event(' in line and 'def button_release_event' not in line:
            indent = len(line) - len(line.lstrip())
            fixed_lines.append(' ' * indent + 'pass  # self.button_release_event(...)')
            if i + 1 < len(lines) and 'guiEvent' in lines[i + 1]:
                skip_next = True
        elif 'self.scroll_event(' in line and 'def scroll_event' not in line:
            indent = len(line) - len(line.lstrip())
            fixed_lines.append(' ' * indent + 'pass  # self.scroll_event(...)')
        else:
            fixed_lines.append(line)
    
    with open(garden_path, 'w') as f:
        f.write('\n'.join(fixed_lines))
    
    print("Garden matplotlib patched successfully!")
PATCH

echo "Setup complete! Run the application with:"
echo "  source .venv/bin/activate"
echo "  python run.py"

