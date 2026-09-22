from __future__ import annotations

import csv
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from .core import predict_formula, predict_csv, OxiMeltError, PredictionResult

class OxiMeltApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("OxiMelt — Oxide Melting Point Calculator")
        self.geometry("1080x760")
        self.minsize(940, 680)

        self.current_result = None
        self.batch_output_path = None

        self._configure_style()
        self._build_ui()

    def _configure_style(self):
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure("Title.TLabel", font=("Segoe UI", 19, "bold"))
        style.configure("Subtitle.TLabel", font=("Segoe UI", 10))
        style.configure("Result.TLabel", font=("Segoe UI", 24, "bold"))
        style.configure("Heading.TLabel", font=("Segoe UI", 11, "bold"))
        style.configure("Status.TLabel", font=("Segoe UI", 10, "bold"))
        style.configure("Treeview", rowheight=26)
        style.configure("Treeview.Heading", font=("Segoe UI", 9, "bold"))

    def _build_ui(self):
        outer = ttk.Frame(self, padding=18)
        outer.pack(fill="both", expand=True)

        ttk.Label(outer, text="OxiMelt", style="Title.TLabel").pack(anchor="w")
        ttk.Label(
            outer,
            text="Formula-only prediction using the frozen global 3P relation with an optional frozen Hf-family local residual transfer.",
            style="Subtitle.TLabel",
        ).pack(anchor="w", pady=(0, 12))

        nb = ttk.Notebook(outer)
        nb.pack(fill="both", expand=True)

        self.single_tab = ttk.Frame(nb, padding=14)
        self.batch_tab = ttk.Frame(nb, padding=14)
        nb.add(self.single_tab, text="Single Formula")
        nb.add(self.batch_tab, text="Batch CSV")

        self._build_single_tab()
        self._build_batch_tab()

        footer = ttk.Label(
            outer,
            text=(
                "Scientific-use software. Predictions are model estimates, not direct experimental measurements. "
                "Check applicability warnings before use."
            ),
            foreground="#555555",
        )
        footer.pack(anchor="w", pady=(10, 0))

    def _build_single_tab(self):
        top = ttk.Frame(self.single_tab)
        top.pack(fill="x")

        ttk.Label(top, text="Chemical formula", style="Heading.TLabel").grid(
            row=0, column=0, sticky="w"
        )
        self.formula_var = tk.StringVar(value="Lu3Al5O12")
        entry = ttk.Entry(top, textvariable=self.formula_var, font=("Consolas", 14), width=36)
        entry.grid(row=1, column=0, sticky="ew", padx=(0, 8), pady=(4, 8))
        entry.bind("<Return>", lambda _e: self._predict_single())

        ttk.Button(top, text="Predict", command=self._predict_single).grid(
            row=1, column=1, padx=4, pady=(4, 8)
        )
        ttk.Button(top, text="Clear", command=self._clear_single).grid(
            row=1, column=2, padx=4, pady=(4, 8)
        )
        ttk.Button(top, text="Export Result", command=self._export_single).grid(
            row=1, column=3, padx=4, pady=(4, 8)
        )
        top.columnconfigure(0, weight=1)

        ttk.Label(
            top,
            text="Examples: Lu3Al5O12   Yb2Ti2O7   Hf0.6Lu0.2Sc0.2O1.8   Lu0.2Yb0.2Hf0.6O1.8",
            foreground="#666666",
        ).grid(row=2, column=0, columnspan=4, sticky="w")

        result_box = ttk.LabelFrame(self.single_tab, text="Prediction", padding=12)
        result_box.pack(fill="x", pady=(14, 10))

        self.result_main = ttk.Label(result_box, text="—", style="Result.TLabel")
        self.result_main.grid(row=0, column=0, sticky="w")
        self.result_c = ttk.Label(result_box, text="", font=("Segoe UI", 13))
        self.result_c.grid(row=1, column=0, sticky="w", pady=(2, 8))

        self.status_label = ttk.Label(result_box, text="Status: —", style="Status.TLabel")
        self.status_label.grid(row=0, column=1, sticky="w", padx=(30, 0))
        self.status_detail = ttk.Label(result_box, text="", wraplength=480, justify="left")
        self.status_detail.grid(row=1, column=1, sticky="nw", padx=(30, 0))

        result_box.columnconfigure(0, weight=1)
        result_box.columnconfigure(1, weight=2)

        self.summary_text = tk.Text(
            self.single_tab, height=8, wrap="word", font=("Consolas", 10)
        )
        self.summary_text.pack(fill="x", pady=(0, 10))
        self.summary_text.configure(state="disabled")

        detail_box = ttk.LabelFrame(self.single_tab, text="Pairwise Calculation Details", padding=8)
        detail_box.pack(fill="both", expand=True)

        columns = (
            "pair", "xi", "xj", "estrain", "echem", "omega", "contribution", "seen"
        )
        self.detail_tree = ttk.Treeview(detail_box, columns=columns, show="headings")
        headings = {
            "pair": "Pair", "xi": "x_i", "xj": "x_j", "estrain": "E_strain",
            "echem": "E_chem", "omega": "Omega_ij (K)",
            "contribution": "x_i x_j Omega_ij (K)", "seen": "Seen in 2-cat training?"
        }
        widths = {
            "pair": 80, "xi": 70, "xj": 70, "estrain": 105, "echem": 105,
            "omega": 120, "contribution": 150, "seen": 150
        }
        for c in columns:
            self.detail_tree.heading(c, text=headings[c])
            self.detail_tree.column(c, width=widths[c], anchor="center")
        ysb = ttk.Scrollbar(detail_box, orient="vertical", command=self.detail_tree.yview)
        self.detail_tree.configure(yscrollcommand=ysb.set)
        self.detail_tree.pack(side="left", fill="both", expand=True)
        ysb.pack(side="right", fill="y")

    def _build_batch_tab(self):
        controls = ttk.Frame(self.batch_tab)
        controls.pack(fill="x")

        ttk.Label(
            controls,
            text="Input CSV must contain a column named: formula",
            style="Heading.TLabel",
        ).pack(side="left")
        ttk.Button(controls, text="Open CSV", command=self._run_batch).pack(side="right", padx=4)
        ttk.Button(
            controls, text="Open Output Folder", command=self._show_batch_output
        ).pack(side="right", padx=4)

        self.batch_status = ttk.Label(self.batch_tab, text="No batch file processed.")
        self.batch_status.pack(fill="x", pady=10)

        cols = ("formula", "global", "corr", "final", "mode", "status", "error")
        self.batch_tree = ttk.Treeview(self.batch_tab, columns=cols, show="headings")
        for c, h, w in [
            ("formula", "Formula", 180),
            ("global", "Global 3P (K)", 120),
            ("corr", "Local correction (K)", 135),
            ("final", "Recommended Tm (K)", 145),
            ("mode", "Prediction mode", 245),
            ("status", "Applicability", 245),
            ("error", "Error", 260),
        ]:
            self.batch_tree.heading(c, text=h)
            self.batch_tree.column(c, width=w, anchor="w")
        self.batch_tree.pack(fill="both", expand=True)

    def _predict_single(self):
        formula = self.formula_var.get().strip()
        try:
            result = predict_formula(formula)
        except Exception as exc:
            self.current_result = None
            self.result_main.config(text="Prediction unavailable")
            self.result_c.config(text="")
            self.status_label.config(text="Status: OUTSIDE CURRENT DOMAIN")
            self.status_detail.config(text=str(exc))
            self._set_summary(f"ERROR\n{exc}")
            self._fill_detail(None)
            return

        self.current_result = result
        self.result_main.config(text=f"{result.predicted_Tm_K:,.1f} K")
        self.result_c.config(text=f"{result.predicted_Tm_C:,.1f} °C  —  Recommended output")
        self.status_label.config(text=f"Status: {result.applicability}")
        self.status_detail.config(text=result.applicability_detail)

        states = ", ".join(
            f"{e}({z:+g})" for e, z in sorted(result.formal_oxidation_states.items())
        )
        fracs = ", ".join(
            f"{e}={x:.4f}" for e, x in sorted(result.cation_fractions.items())
        )
        warning_text = "\n".join(f"NOTE: {w}" for w in result.warnings) or "Notes: none"
        local_lines = ""
        if result.hf_local_domain_eligible:
            local_lines += (
                f"Hf-family support coverage: {result.hf_local_support_coverage:.1%}\n"
                f"Local correction applied: {'Yes' if result.hf_local_correction_applied else 'No'}\n"
            )
            for d in result.local_support_details:
                sx = (
                    f"x_Hf={d.support_x_Hf:.4f}" if d.support_x_Hf is not None
                    else (f"x_Hf=[{d.support_x_Hf_low:.4f}, {d.support_x_Hf_high:.4f}]"
                          if d.support_x_Hf_low is not None else "")
                )
                local_lines += (
                    f"  {d.dopant}: share={d.dopant_share:.3f}, {d.mode}, "
                    f"residual={d.residual_K:+.3f} K, {sx}\n"
                )
        summary = (
            f"Formula: {result.formula}\n"
            f"Prediction mode: {result.prediction_mode}\n"
            f"Formal oxidation states: {states}\n"
            f"Cation fractions: {fracs}\n"
            f"T0 (simple-oxide reference): {result.T0_K:.3f} K\n"
            f"Global pair correction: {result.pair_correction_K:+.3f} K\n"
            f"Global 3P prediction: {result.global_Tm_K:.3f} K\n"
            f"Family-specific local correction: {result.local_correction_K:+.3f} K\n"
            f"Recommended predicted Tm: {result.predicted_Tm_K:.3f} K "
            f"({result.predicted_Tm_C:.3f} °C)\n"
            f"{local_lines}"
            f"{warning_text}"
        )
        self._set_summary(summary)
        self._fill_detail(result)

    def _set_summary(self, text):
        self.summary_text.configure(state="normal")
        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("1.0", text)
        self.summary_text.configure(state="disabled")

    def _fill_detail(self, result):
        for item in self.detail_tree.get_children():
            self.detail_tree.delete(item)
        if result is None:
            return
        for p in result.pair_contributions:
            self.detail_tree.insert(
                "", "end",
                values=(
                    p.pair, f"{p.x_i:.5f}", f"{p.x_j:.5f}",
                    f"{p.E_strain:.6f}", f"{p.E_chem:.6f}",
                    f"{p.Omega_K:.3f}", f"{p.weighted_contribution_K:+.3f}",
                    "Yes" if p.seen_in_two_cation_training else "No"
                )
            )

    def _clear_single(self):
        self.formula_var.set("")
        self.current_result = None
        self.result_main.config(text="—")
        self.result_c.config(text="")
        self.status_label.config(text="Status: —")
        self.status_detail.config(text="")
        self._set_summary("")
        self._fill_detail(None)

    def _export_single(self):
        if self.current_result is None:
            messagebox.showinfo("OxiMelt", "Run a prediction first.")
            return
        path = filedialog.asksaveasfilename(
            title="Export result",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile=f"{self.current_result.formula}_oximelt_prediction.csv",
        )
        if not path:
            return

        r = self.current_result
        fields = list(r.to_flat_dict().keys())
        with open(path, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerow(r.to_flat_dict())
        messagebox.showinfo("OxiMelt", f"Result exported to:\n{path}")

    def _run_batch(self):
        path = filedialog.askopenfilename(
            title="Select batch CSV",
            filetypes=[("CSV files", "*.csv")]
        )
        if not path:
            return
        out = str(Path(path).with_name(Path(path).stem + "_oximelt_predictions.csv"))
        try:
            n = predict_csv(path, out)
        except Exception as exc:
            messagebox.showerror("Batch prediction failed", str(exc))
            return
        self.batch_output_path = Path(out)
        self.batch_status.config(text=f"Processed {n} formula(s). Output: {out}")
        self._load_batch_output(out)

    def _load_batch_output(self, path):
        for item in self.batch_tree.get_children():
            self.batch_tree.delete(item)
        with open(path, newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                self.batch_tree.insert(
                    "", "end",
                    values=(
                        row["formula"],
                        row.get("global_Tm_K", ""),
                        row.get("local_correction_K", ""),
                        row.get("predicted_Tm_K", ""),
                        row.get("prediction_mode", ""),
                        row.get("applicability", ""),
                        row.get("error", ""),
                    )
                )

    def _show_batch_output(self):
        if not self.batch_output_path:
            messagebox.showinfo("OxiMelt", "No batch output is available yet.")
            return
        messagebox.showinfo(
            "Batch output",
            f"Output file:\n{self.batch_output_path}\n\n"
            "Open this path with your system file manager."
        )

def main():
    app = OxiMeltApp()
    app.mainloop()

if __name__ == "__main__":
    main()
