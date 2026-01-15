import tkinter as tk
from tkinter import ttk
from typing import Callable, List, Optional

from point_mass import PointMass
from si_prefix import prefix_labels, prefix_multiplier
from torsional_spring import TorsionalSpring
from translational_spring import TranslationalSpring


class _ListItem:
    def __init__(self, title: str, on_select: Callable[[Optional[int]], None]):
        self.title = title
        self.on_select = on_select
        self.items: List = []
        self.listbox: Optional[tk.Listbox] = None

    def build(self, master: tk.Widget):
        frame = ttk.Frame(master)
        ttk.Label(frame, text=self.title).pack(anchor="w")
        self.listbox = tk.Listbox(frame, height=4, exportselection=False)
        self.listbox.pack(fill="both", expand=True)
        self.listbox.bind("<<ListboxSelect>>", self._on_select)
        return frame

    def _on_select(self, _event=None):
        if not self.listbox:
            return
        sel = self.listbox.curselection()
        idx = sel[0] if sel else None
        self.on_select(idx)

    def refresh(self, fmt: Callable[[int, object], str]):
        if not self.listbox:
            return
        self.listbox.delete(0, tk.END)
        for i, item in enumerate(self.items):
            self.listbox.insert(tk.END, fmt(i, item))


class AttachmentsFrame(ttk.LabelFrame):
    def __init__(self, master: tk.Widget, length_getter: Callable[[], float], on_change: Callable[[], None]):
        super().__init__(master, text="Attachments")
        self.length_getter = length_getter
        self.on_change = on_change

        self.masses = _ListItem("Point masses", self._select_mass)
        self.trans_springs = _ListItem("Trans springs", self._select_trans)
        self.tors_springs = _ListItem("Tors springs", self._select_tors)

        list_frame = ttk.Frame(self)
        self.mass_widget = self.masses.build(list_frame)
        self.trans_widget = self.trans_springs.build(list_frame)
        self.tors_widget = self.tors_springs.build(list_frame)
        self.mass_widget.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)
        self.trans_widget.grid(row=0, column=1, sticky="nsew", padx=4, pady=4)
        self.tors_widget.grid(row=0, column=2, sticky="nsew", padx=4, pady=4)
        list_frame.columnconfigure((0, 1, 2), weight=1)
        list_frame.grid(row=0, column=0, columnspan=2, sticky="nsew")

        form = ttk.Frame(self)
        ttk.Label(form, text="Position [m]").grid(row=0, column=0, sticky="w")
        self.pos_scale = tk.Scale(
            form, from_=0.0, to=1.0, orient="horizontal", resolution=0.001, command=self._pos_changed
        )
        self.pos_scale.grid(row=0, column=1, sticky="ew", padx=4, pady=2)

        ttk.Label(form, text="Mass [kg]").grid(row=1, column=0, sticky="w")
        self.mass_var = tk.StringVar(value="1.0")
        ttk.Entry(form, textvariable=self.mass_var, width=10).grid(row=1, column=1, sticky="ew", padx=4, pady=2)
        self.mass_prefix_var = tk.StringVar(value="kg")
        ttk.OptionMenu(
            form,
            self.mass_prefix_var,
            self.mass_prefix_var.get(),
            *prefix_labels("kg"),
        ).grid(row=1, column=2, sticky="w", padx=4, pady=2)

        ttk.Label(form, text="J [kg*m^2]").grid(row=2, column=0, sticky="w")
        self.j_var = tk.StringVar(value="0.0")
        ttk.Entry(form, textvariable=self.j_var, width=10).grid(row=2, column=1, sticky="ew", padx=4, pady=2)
        self.j_prefix_var = tk.StringVar(value="kg*m^2")
        ttk.OptionMenu(
            form,
            self.j_prefix_var,
            self.j_prefix_var.get(),
            *prefix_labels("kg*m^2"),
        ).grid(row=2, column=2, sticky="w", padx=4, pady=2)

        ttk.Label(form, text="k trans [N/m]").grid(row=3, column=0, sticky="w")
        self.kt_var = tk.StringVar(value="1e5")
        ttk.Entry(form, textvariable=self.kt_var, width=10).grid(row=3, column=1, sticky="ew", padx=4, pady=2)
        self.kt_prefix_var = tk.StringVar(value="N/m")
        ttk.OptionMenu(
            form,
            self.kt_prefix_var,
            self.kt_prefix_var.get(),
            *prefix_labels("N/m"),
        ).grid(row=3, column=2, sticky="w", padx=4, pady=2)

        ttk.Label(form, text="k tors [N*m/rad]").grid(row=4, column=0, sticky="w")
        self.kq_var = tk.StringVar(value="1e4")
        ttk.Entry(form, textvariable=self.kq_var, width=10).grid(row=4, column=1, sticky="ew", padx=4, pady=2)
        self.kq_prefix_var = tk.StringVar(value="N*m/rad")
        ttk.OptionMenu(
            form,
            self.kq_prefix_var,
            self.kq_prefix_var.get(),
            *prefix_labels("N*m/rad"),
        ).grid(row=4, column=2, sticky="w", padx=4, pady=2)

        btns = ttk.Frame(form)
        ttk.Button(btns, text="Add mass", command=self._add_mass).grid(row=0, column=0, padx=2)
        ttk.Button(btns, text="Del mass", command=self._del_mass).grid(row=0, column=1, padx=2)
        ttk.Button(btns, text="Add trans spring", command=self._add_trans).grid(row=1, column=0, padx=2)
        ttk.Button(btns, text="Del trans spring", command=self._del_trans).grid(row=1, column=1, padx=2)
        ttk.Button(btns, text="Add tors spring", command=self._add_tors).grid(row=2, column=0, padx=2)
        ttk.Button(btns, text="Del tors spring", command=self._del_tors).grid(row=2, column=1, padx=2)
        btns.grid(row=5, column=0, columnspan=2, pady=4)

        form.columnconfigure(1, weight=1)
        form.grid(row=1, column=0, sticky="nsew", padx=4, pady=4)

        self.selected_mass: Optional[int] = None
        self.selected_trans: Optional[int] = None
        self.selected_tors: Optional[int] = None

    def _pos_changed(self, _value=None):
        pos = float(self.pos_scale.get()) * self.length_getter()
        if self.selected_mass is not None:
            pm = self.masses.items[self.selected_mass]
            self.masses.items[self.selected_mass] = PointMass(position=pos, mass=pm.mass, rotary_inertia=pm.rotary_inertia)
        if self.selected_trans is not None:
            sp = self.trans_springs.items[self.selected_trans]
            self.trans_springs.items[self.selected_trans] = TranslationalSpring(position=pos, k=sp.k)
        if self.selected_tors is not None:
            sp = self.tors_springs.items[self.selected_tors]
            self.tors_springs.items[self.selected_tors] = TorsionalSpring(position=pos, k=sp.k)
        self._refresh_lists()
        self.on_change()

    def _add_mass(self):
        try:
            mass = float(self.mass_var.get()) * prefix_multiplier(self.mass_prefix_var.get(), "kg")
            j = float(self.j_var.get()) * prefix_multiplier(self.j_prefix_var.get(), "kg*m^2")
        except ValueError:
            return
        pos = float(self.pos_scale.get()) * self.length_getter()
        self.masses.items.append(PointMass(position=pos, mass=mass, rotary_inertia=j))
        self._refresh_lists()
        self.on_change()

    def _del_mass(self):
        if self.selected_mass is None:
            return
        self.masses.items.pop(self.selected_mass)
        self.selected_mass = None
        self._refresh_lists()
        self.on_change()

    def _add_trans(self):
        try:
            k = float(self.kt_var.get()) * prefix_multiplier(self.kt_prefix_var.get(), "N/m")
        except ValueError:
            return
        pos = float(self.pos_scale.get()) * self.length_getter()
        self.trans_springs.items.append(TranslationalSpring(position=pos, k=k))
        self._refresh_lists()
        self.on_change()

    def _del_trans(self):
        if self.selected_trans is None:
            return
        self.trans_springs.items.pop(self.selected_trans)
        self.selected_trans = None
        self._refresh_lists()
        self.on_change()

    def _add_tors(self):
        try:
            k = float(self.kq_var.get()) * prefix_multiplier(self.kq_prefix_var.get(), "N*m/rad")
        except ValueError:
            return
        pos = float(self.pos_scale.get()) * self.length_getter()
        self.tors_springs.items.append(TorsionalSpring(position=pos, k=k))
        self._refresh_lists()
        self.on_change()

    def _del_tors(self):
        if self.selected_tors is None:
            return
        self.tors_springs.items.pop(self.selected_tors)
        self.selected_tors = None
        self._refresh_lists()
        self.on_change()

    def _select_mass(self, idx: Optional[int]):
        self.selected_mass = idx
        if idx is not None:
            pm = self.masses.items[idx]
            self.mass_var.set(str(pm.mass))
            self.j_var.set(str(pm.rotary_inertia))
            self.mass_prefix_var.set("kg")
            self.j_prefix_var.set("kg*m^2")
            self.pos_scale.set(pm.position / max(self.length_getter(), 1e-9))

    def _select_trans(self, idx: Optional[int]):
        self.selected_trans = idx
        if idx is not None:
            sp = self.trans_springs.items[idx]
            self.kt_var.set(str(sp.k))
            self.kt_prefix_var.set("N/m")
            self.pos_scale.set(sp.position / max(self.length_getter(), 1e-9))

    def _select_tors(self, idx: Optional[int]):
        self.selected_tors = idx
        if idx is not None:
            sp = self.tors_springs.items[idx]
            self.kq_var.set(str(sp.k))
            self.kq_prefix_var.set("N*m/rad")
            self.pos_scale.set(sp.position / max(self.length_getter(), 1e-9))

    def _refresh_lists(self):
        self.masses.refresh(lambda i, m: f"{i}: x={m.position:.2f} m m={m.mass:g} J={m.rotary_inertia:g}")
        self.trans_springs.refresh(lambda i, k: f"{i}: x={k.position:.2f} m k={k.k:g}")
        self.tors_springs.refresh(lambda i, k: f"{i}: x={k.position:.2f} m k={k.k:g}")

    def get_items(self):
        return (
            list(self.masses.items),
            list(self.trans_springs.items),
            list(self.tors_springs.items),
        )
