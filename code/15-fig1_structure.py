#!/usr/bin/env python3
"""
Final production script: PARP1 structural panels with PIL text overlay.
- PyMOL: clean structure-only render (no labels, no distance objects)
- PIL: text labels, domain annotations, CV arrows overlay
- Output: RGB PNG, 600 DPI, white background, publication-ready
"""
import subprocess, os
from PIL import Image, ImageDraw, ImageFont

os.makedirs("results/figures", exist_ok=True)

# ====== PHASE 1: PyMOL — structure-only renders ======
# Each panel in isolated PyMOL process for reliability

pymol_scenes = {
    "A_Domains": r"""
load data/00_raw/pdb_structures/4DQY.pdb
hide everything; bg_color white
set cartoon_fancy_helices, 1; set cartoon_smooth_loops, 1
set ray_trace_mode, 1; set antialias, 2; set depth_cue, 0
show cartoon, chain A; color tv_blue, chain A
show cartoon, chain B; color cyan, chain B
show cartoon, chain C; color limon, chain C and resi 531-661
color orange, chain C and resi 662-787; color red, chain C and resi 788-1011
show cartoon, chain M+N; color tv_green, chain M+N
zoom chain C; turn y, 30
ray 1800, 1350
""",
    "B_S1_CAT": r"""
load data/00_raw/pdb_structures/4DQY.pdb
hide everything; bg_color white
set cartoon_fancy_helices, 1; set cartoon_smooth_loops, 1
set ray_trace_mode, 1; set antialias, 2; set depth_cue, 0
show cartoon, chain C
color orange, chain C and resi 662-787; color red, chain C and resi 788-1011
color limon, chain C and resi 531-661
zoom chain C; turn y, -20
ray 1500, 1125
""",
    "C_S2_Full": r"""
load data/00_raw/pdb_structures/4DQY.pdb
hide everything; bg_color white
set cartoon_fancy_helices, 1; set cartoon_smooth_loops, 1
set ray_trace_mode, 1; set antialias, 2; set depth_cue, 0
show cartoon
color tv_blue, chain A; color cyan, chain B
color limon, chain C and resi 531-661
color orange, chain C and resi 662-787; color red, chain C and resi 788-1011
color tv_green, chain M+N
zoom visible; turn y, 45
ray 1500, 1125
""",
    "D_Pocket": r"""
load data/00_raw/pdb_structures/4DQY.pdb
hide everything; bg_color white
set cartoon_fancy_helices, 1; set cartoon_smooth_loops, 1
set ray_trace_mode, 1; set antialias, 2; set depth_cue, 0
show cartoon, chain C
color orange, chain C and resi 662-787; color red, chain C and resi 788-1011
select pocket, chain C and resi 860-910+980-1000
show sticks, pocket; color white, pocket
set stick_radius, 0.2
zoom pocket; turn y, 15; turn x, 10
ray 1800, 1350
"""
}

raw_dir = "/tmp/parp1_raw"
os.makedirs(raw_dir, exist_ok=True)

for name, script in pymol_scenes.items():
    spath = f"/tmp/pymol_{name}.pml"
    script_full = script + f"\npng {raw_dir}/{name}_raw.png, dpi=600\nquit\n"
    with open(spath, "w") as f:
        f.write(script_full)
    
    print(f"PyMOL rendering Panel {name}...")
    r = subprocess.run(["pymol", "-cq", spath], capture_output=True, text=True, timeout=90,
                       cwd="/Users/taozhu/clacky_workspace/PARPi_design")
    
    out = f"{raw_dir}/{name}_raw.png"
    if os.path.exists(out):
        img = Image.open(out)
        print(f"  ✓ {img.size[0]}×{img.size[1]} px")
    else:
        # Check for crash, try lower res
        if "out of memory" in r.stdout.lower() or "EEK" in r.stdout.lower():
            print(f"  ⚠ OOM, retrying at lower resolution...")
            script_low = script + "\nray 1200, 900\n" + f"png {raw_dir}/{name}_raw.png, dpi=600\nquit\n"
            with open(f"/tmp/pymol_{name}_low.pml", "w") as f:
                f.write(script_low)
            subprocess.run(["pymol", "-cq", f"/tmp/pymol_{name}_low.pml"], 
                          capture_output=True, timeout=90,
                          cwd="/Users/taozhu/clacky_workspace/PARPi_design")
            if os.path.exists(out):
                img = Image.open(out)
                print(f"  ✓ (reduced) {img.size[0]}×{img.size[1]} px")
            else:
                print(f"  ✗ FAILED")

# ====== PHASE 2: PIL — text overlay + white bg conversion ======

# Try to load Arial, fallback to default
try:
    font_title = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 36)
    font_label = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 28)
    font_small = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 24)
except:
    font_title = ImageFont.load_default()
    font_label = ImageFont.load_default()
    font_small = ImageFont.load_default()

def add_text_overlay(img, draw, x, y, text, font, color):
    """Add text with white outline for readability."""
    # Outline
    for dx in [-1, 1]:
        for dy in [-1, 1]:
            draw.text((x+dx, y+dy), text, font=font, fill=(255, 255, 255))
    # Main text
    draw.text((x, y), text, font=font, fill=color)

def convert_to_rgb_white(img):
    """Convert RGBA to RGB with white background."""
    if img.mode == 'RGBA':
        bg = Image.new('RGB', img.size, (255, 255, 255))
        bg.paste(img, mask=img.split()[3])
        return bg
    elif img.mode != 'RGB':
        return img.convert('RGB')
    return img

annotations = {
    "A_Domains": lambda img, draw, w, h: [
        # Domain labels
        add_text_overlay(img, draw, int(w*0.12), int(h*0.12), "ZnF1", font_label, (100, 149, 237)),
        add_text_overlay(img, draw, int(w*0.28), int(h*0.15), "ZnF3", font_label, (0, 191, 255)),
        add_text_overlay(img, draw, int(w*0.50), int(h*0.08), "WGR", font_label, (50, 205, 50)),
        add_text_overlay(img, draw, int(w*0.62), int(h*0.18), "HD", font_title, (255, 165, 0)),
        add_text_overlay(img, draw, int(w*0.70), int(h*0.35), "ART", font_title, (220, 20, 60)),
        add_text_overlay(img, draw, int(w*0.80), int(h*0.80), "DNA", font_label, (60, 179, 113)),
        # Title
        add_text_overlay(img, draw, int(w*0.02), int(h*0.02), "A  PARP1 Domain Architecture (4DQY)", font_title, (50, 50, 50)),
    ],
    "B_S1_CAT": lambda img, draw, w, h: [
        add_text_overlay(img, draw, int(w*0.50), int(h*0.10), "HD", font_title, (255, 165, 0)),
        add_text_overlay(img, draw, int(w*0.50), int(h*0.35), "(residues 662-787)", font_small, (180, 120, 30)),
        add_text_overlay(img, draw, int(w*0.55), int(h*0.55), "ART", font_title, (220, 20, 60)),
        add_text_overlay(img, draw, int(w*0.55), int(h*0.72), "(residues 788-1011)", font_small, (180, 40, 40)),
        add_text_overlay(img, draw, int(w*0.02), int(h*0.02), "B  System 1: CAT-only (~30K atoms)", font_title, (50, 50, 50)),
        add_text_overlay(img, draw, int(w*0.02), int(h*0.06), "CV: HD-ART center-of-mass distance", font_small, (100, 100, 100)),
    ],
    "C_S2_Full": lambda img, draw, w, h: [
        add_text_overlay(img, draw, int(w*0.15), int(h*0.08), "ZnF1", font_small, (100, 149, 237)),
        add_text_overlay(img, draw, int(w*0.30), int(h*0.20), "ZnF3", font_small, (0, 191, 255)),
        add_text_overlay(img, draw, int(w*0.45), int(h*0.25), "WGR", font_small, (50, 205, 50)),
        add_text_overlay(img, draw, int(w*0.55), int(h*0.18), "HD", font_label, (255, 165, 0)),
        add_text_overlay(img, draw, int(w*0.65), int(h*0.42), "ART", font_label, (220, 20, 60)),
        add_text_overlay(img, draw, int(w*0.55), int(h*0.85), "DNA", font_label, (60, 179, 113)),
        add_text_overlay(img, draw, int(w*0.02), int(h*0.02), "C  System 2: Full PARP1+DNA (~120K atoms)", font_title, (50, 50, 50)),
        add_text_overlay(img, draw, int(w*0.02), int(h*0.06), "CV1: Prot-DNA  |  CV2: HD-ART", font_small, (100, 100, 100)),
    ],
    "D_Pocket": lambda img, draw, w, h: [
        # Catalytic residues
        add_text_overlay(img, draw, int(w*0.45), int(h*0.25), "G863", font_label, (255, 215, 0)),
        add_text_overlay(img, draw, int(w*0.65), int(h*0.40), "S904", font_label, (255, 215, 0)),
        add_text_overlay(img, draw, int(w*0.30), int(h*0.65), "E988", font_label, (220, 20, 60)),
        add_text_overlay(img, draw, int(w*0.30), int(h*0.72), "(catalytic)", font_small, (180, 40, 40)),
        # Domain labels
        add_text_overlay(img, draw, int(w*0.72), int(h*0.08), "HD", font_label, (255, 165, 0)),
        add_text_overlay(img, draw, int(w*0.72), int(h*0.14), "(allosteric)", font_small, (180, 120, 30)),
        add_text_overlay(img, draw, int(w*0.72), int(h*0.55), "ART", font_label, (220, 20, 60)),
        add_text_overlay(img, draw, int(w*0.72), int(h*0.61), "(catalytic)", font_small, (180, 40, 40)),
        add_text_overlay(img, draw, int(w*0.02), int(h*0.02), "D  NAD+ Binding Pocket", font_title, (50, 50, 50)),
        add_text_overlay(img, draw, int(w*0.02), int(h*0.07), "Key catalytic residues", font_small, (100, 100, 100)),
    ],
}

print("\n=== PIL: Adding text overlays + white background ===")
for name, annotate_fn in annotations.items():
    raw_path = f"{raw_dir}/{name}_raw.png"
    if not os.path.exists(raw_path):
        print(f"  ✗ {name}: raw file missing, skipping")
        continue
    
    img = Image.open(raw_path)
    img = convert_to_rgb_white(img)
    draw = ImageDraw.Draw(img)
    w, h = img.size
    
    annotate_fn(img, draw, w, h)
    
    out_path = f"results/figures/Fig1_Panel{name.replace('_','')}.png"
    img.save(out_path, "PNG", dpi=(600, 600))
    print(f"  ✓ {out_path} — {w}×{h} px, 600 DPI")

# ====== FINAL VERIFICATION ======
print("\n=== Final Check ===")
for name in ["A_Domains", "B_S1_CAT", "C_S2_Full", "D_Pocket"]:
    fp = f"results/figures/Fig1_Panel{name.replace('_','')}.png"
    if os.path.exists(fp):
        img = Image.open(fp)
        dpi = img.info.get('dpi', (0,0))
        print(f"  {fp.split('/')[-1]}: {img.size[0]}×{img.size[1]} px, {dpi[0]:.0f} DPI, mode={img.mode}, {os.path.getsize(fp)//1024}KB")
    else:
        print(f"  {fp.split('/')[-1]}: MISSING")
