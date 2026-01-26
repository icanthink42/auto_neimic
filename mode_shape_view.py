import tkinter as tk
from tkinter import ttk
from typing import Optional, Sequence

import numpy as np


class ModeShapeView(ttk.LabelFrame):
    def __init__(self, master: tk.Widget):
        super().__init__(master, text="Mode Shapes")
        self.canvas = tk.Canvas(
            self, width=800, height=300, bg="white", highlightthickness=1, highlightbackground="#ccc"
        )
        self.canvas.pack(fill="both", expand=True, padx=4, pady=4)
        
        # Control panel
        controls = ttk.Frame(self)
        controls.pack(fill="x", padx=4, pady=2)
        
        ttk.Label(controls, text="Mode:").pack(side="left", padx=2)
        self.mode_var = tk.IntVar(value=1)
        self.mode_spinbox = ttk.Spinbox(
            controls, from_=1, to=1, textvariable=self.mode_var,
            width=5, command=self._on_mode_change
        )
        self.mode_spinbox.pack(side="left", padx=2)
        
        self.type_var = tk.StringVar(value="bending")
        ttk.Radiobutton(
            controls, text="Bending", variable=self.type_var, value="bending",
            command=self._on_type_change
        ).pack(side="left", padx=4)
        ttk.Radiobutton(
            controls, text="Torsion", variable=self.type_var, value="torsion",
            command=self._on_type_change
        ).pack(side="left", padx=4)
        
        self.freq_label = ttk.Label(controls, text="Frequency: —")
        self.freq_label.pack(side="left", padx=8)
        
        # Data storage
        self.beam_length = 2.0
        self.n_nodes = 21
        self.bend_freqs = np.array([])
        self.tors_freqs = np.array([])
        self.bend_modes = np.zeros((0, 0))
        self.tors_modes = np.zeros((0, 0))
        self.bending_fixed = []
        self.torsion_fixed = []
        
        self.margin = 40
        
        self.canvas.bind("<Configure>", lambda _e: self._draw())
    
    def update_modes(
        self,
        beam_length: float,
        n_nodes: int,
        bend_freqs: np.ndarray,
        tors_freqs: np.ndarray,
        bend_modes: np.ndarray,
        tors_modes: np.ndarray,
        bending_fixed: Sequence[int],
        torsion_fixed: Sequence[int],
    ):
        """Update the mode shape data and redraw."""
        self.beam_length = beam_length
        self.n_nodes = n_nodes
        self.bend_freqs = bend_freqs
        self.tors_freqs = tors_freqs
        self.bend_modes = bend_modes
        self.tors_modes = tors_modes
        self.bending_fixed = list(bending_fixed)
        self.torsion_fixed = list(torsion_fixed)
        
        # Update spinbox range
        n_bend = len(bend_freqs)
        n_tors = len(tors_freqs)
        max_modes = max(n_bend, n_tors, 1)
        self.mode_spinbox.config(to=max_modes)
        
        # Reset to mode 1 if current mode is out of range
        if self.mode_var.get() > max_modes:
            self.mode_var.set(1)
        
        self._draw()
    
    def _on_mode_change(self):
        self._draw()
    
    def _on_type_change(self):
        self._draw()
    
    def _dims(self):
        w = max(300, self.canvas.winfo_width())
        h = max(300, self.canvas.winfo_height())
        return w, h
    
    def _draw(self):
        self.canvas.delete("all")
        
        mode_type = self.type_var.get()
        mode_idx = self.mode_var.get() - 1  # Convert to 0-based index
        
        if mode_type == "bending":
            freqs = self.bend_freqs
            modes = self.bend_modes
            fixed_dofs = self.bending_fixed
        else:
            freqs = self.tors_freqs
            modes = self.tors_modes
            fixed_dofs = self.torsion_fixed
        
        if len(freqs) == 0 or mode_idx >= len(freqs):
            self.freq_label.config(text="Frequency: —")
            w, h = self._dims()
            self.canvas.create_text(
                w // 2, h // 2,
                text="No mode shapes available\n(Run analysis with fixed boundary conditions)",
                fill="#999", font=("Arial", 12)
            )
            return
        
        # Update frequency label
        freq = freqs[mode_idx]
        self.freq_label.config(text=f"Frequency: {freq:.3f} Hz")
        
        # Get the mode shape (reduced DOFs)
        mode_reduced = modes[:, mode_idx]
        
        # Reconstruct full DOF vector
        if mode_type == "bending":
            n_dof = 2 * self.n_nodes
            mode_full = np.zeros(n_dof)
            free_idx = [i for i in range(n_dof) if i not in fixed_dofs]
            mode_full[free_idx] = mode_reduced
            
            # Extract displacement (every other DOF starting at 0)
            displacement = mode_full[::2]
        else:
            # Torsion: one DOF per node
            mode_full = np.zeros(self.n_nodes)
            free_idx = [i for i in range(self.n_nodes) if i not in fixed_dofs]
            mode_full[free_idx] = mode_reduced
            displacement = mode_full
        
        self._draw_mode_shape(displacement, mode_type)
    
    def _draw_mode_shape(self, displacement: np.ndarray, mode_type: str):
        """Draw the mode shape."""
        w, h = self._dims()
        m = self.margin
        x0, x1 = m, w - m
        y_center = h // 2
        
        # Scale displacement for visualization
        max_disp = np.max(np.abs(displacement))
        if max_disp < 1e-10:
            max_disp = 1.0
        
        scale_factor = (h - 2 * m) * 0.35 / max_disp
        
        # X positions along beam
        x_beam = np.linspace(0, self.beam_length, len(displacement))
        
        # Draw reference (undeformed) beam
        self.canvas.create_line(
            x0, y_center, x1, y_center,
            fill="#ccc", width=2, dash=(4, 4)
        )
        
        # Draw deformed shape
        points = []
        for i, (x_pos, disp) in enumerate(zip(x_beam, displacement)):
            cx = x0 + (x_pos / self.beam_length) * (x1 - x0)
            
            if mode_type == "bending":
                cy = y_center - disp * scale_factor
            else:
                # For torsion, show as rotation (simplified as vertical displacement)
                cy = y_center - disp * scale_factor
            
            points.append((cx, cy))
        
        # Draw mode shape line
        if len(points) > 1:
            flat_points = [coord for point in points for coord in point]
            self.canvas.create_line(
                *flat_points,
                fill="#e74c3c", width=3, smooth=True
            )
            
            # Draw nodes
            for cx, cy in points:
                self.canvas.create_oval(
                    cx - 3, cy - 3, cx + 3, cy + 3,
                    fill="#c0392b", outline=""
                )
        
        # Draw axis labels
        self.canvas.create_text(
            x0, y_center + 20, text="0", anchor="w", fill="#666"
        )
        self.canvas.create_text(
            x1, y_center + 20, text=f"{self.beam_length:g} m", anchor="e", fill="#666"
        )
        
        # Draw title
        mode_num = self.mode_var.get()
        title = f"Mode {mode_num} - {mode_type.capitalize()}"
        self.canvas.create_text(
            w // 2, 15, text=title, fill="#333", font=("Arial", 12, "bold")
        )

