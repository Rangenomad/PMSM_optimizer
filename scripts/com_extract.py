"""COM-based data extraction for AEDT — version-agnostic, stable since AEDT 2015.

Use this when pyaedt's post-processing API has gRPC version mismatches with
the installed AEDT. COM interface is AEDT's native protocol and doesn't change.

Usage:
    from scripts.com_extract import extract_transient
    time_s, data_dict = extract_transient(
        project_path='pmsm_projects/2026-07-27_bemf_3000rpm_final',
        design_name='5_Partial_motor_TR',
        expressions=['InducedVoltage(Phase_A)', 'InducedVoltage(Phase_B)', 'InducedVoltage(Phase_C)'],
    )
"""

import sys
import csv
from pathlib import Path

import numpy as np
import win32com.client


def extract_transient(project_path, design_name='5_Partial_motor_TR',
                      expressions=None, setup_name=None):
    """Extract transient solution data via COM (version-agnostic).

    Parameters
    ----------
    project_path : str or Path
        Path to the project directory (containing Prius_2D_Practice.aedt).
    design_name : str
        AEDT design name (e.g. '5_Partial_motor_TR').
    expressions : list of str, optional
        Expressions to extract. Defaults to 3-phase induced voltages.
    setup_name : str, optional
        Setup name. Auto-detected if None.

    Returns
    -------
    time_s : np.ndarray
        Time vector in seconds.
    data : dict
        {expression_name: np.ndarray} of extracted waveforms.
    """
    if expressions is None:
        expressions = [
            'InducedVoltage(Phase_A)',
            'InducedVoltage(Phase_B)',
            'InducedVoltage(Phase_C)',
        ]

    project_dir = Path(project_path)
    aedt_file = project_dir / 'Prius_2D_Practice.aedt'
    if not aedt_file.exists():
        # Try direct path
        aedt_file = Path(project_path)
    project_file = str(aedt_file.resolve())

    # 1. Connect via COM
    oApp = win32com.client.Dispatch('Ansoft.ElectronicsDesktop')
    desktop = oApp.GetAppDesktop()

    # 2. Open project (skip if already open)
    project_name = aedt_file.stem  # 'Prius_2D_Practice'
    try:
        existing = desktop.GetActiveProject()
        if existing is not None and existing.GetName() == project_name:
            pass  # Already open
        else:
            desktop.OpenProject(project_file)
    except Exception:
        try:
            desktop.OpenProject(project_file)
        except Exception as e2:
            # If still locked, try to find it among open projects
            print(f'[COM] OpenProject failed ({e2}), searching open projects...')
            for proj_name in desktop.GetProjects():
                if proj_name == project_name:
                    desktop.SetActiveProject(project_name)
                    break
            else:
                raise RuntimeError(f'Cannot open or find project: {project_file}')

    proj = desktop.GetActiveProject()
    design = proj.GetDesign(design_name)
    reporter = design.GetModule("ReportSetup")

    # 3. Discover setup name
    if setup_name is None:
        try:
            solutions = reporter.GetAvailableSolutions("Transient")
        except Exception:
            # Try without argument (older AEDT versions)
            solutions = reporter.GetAvailableSolutions()
        if not solutions:
            # Fallback to common default
            setup_name = "Setup1 : Transient"
            print(f'[COM] Using default setup: {setup_name}')
        else:
            setup_name = solutions[0]
            print(f'[COM] Auto-detected setup: {setup_name}')

    # 4. Extract data via GetSolutionDataPerVariation (direct, no report needed)
    result = reporter.GetSolutionDataPerVariation(
        "Transient",        # reportCategory
        setup_name,         # setupName
        [],                 # sweepList
        [],                 # variationList
        list(expressions),  # expressions
    )

    # 5. Parse the COM result object
    # GetSolutionDataPerVariation returns a COM object with:
    #   .GetSweepValues() -> list of sweep values (Time in seconds)
    #   .GetRealDataValues(expression) -> real part values
    #   .GetImagDataValues(expression) -> imaginary part values

    time_vals = np.array(list(result.GetSweepValues()), dtype=float)
    data = {}
    for expr in expressions:
        real_vals = np.array(list(result.GetRealDataValues(expr)), dtype=float)
        data[expr] = real_vals

    # 6. Cleanup
    desktop.CloseProject(project_name)

    return time_vals, data


def extract_magnetostatic(project_path, design_name='4_Partial_motor_MS2',
                          expressions=None, setup_name=None):
    """Extract magnetostatic solution data via COM.

    Parameters are the same as extract_transient().
    """
    if expressions is None:
        expressions = [
            'FluxLinkage(Phase_A)', 'FluxLinkage(Phase_B)', 'FluxLinkage(Phase_C)',
        ]
    # Same flow but with "Magnetostatic" category
    project_dir = Path(project_path)
    aedt_file = project_dir / 'Prius_2D_Practice.aedt'
    if not aedt_file.exists():
        aedt_file = Path(project_path)
    project_file = str(aedt_file.resolve())
    project_name = aedt_file.stem

    oApp = win32com.client.Dispatch('Ansoft.ElectronicsDesktop')
    desktop = oApp.GetAppDesktop()

    try:
        existing = desktop.GetActiveProject()
        if existing is None or existing.GetName() != project_name:
            desktop.OpenProject(project_file)
    except Exception:
        desktop.OpenProject(project_file)

    proj = desktop.GetActiveProject()
    design = proj.GetDesign(design_name)
    reporter = design.GetModule("ReportSetup")

    if setup_name is None:
        solutions = reporter.GetAvailableSolutions()
        if not solutions:
            raise RuntimeError(f'No solutions found in design {design_name}')
        setup_name = solutions[0]

    result = reporter.GetSolutionDataPerVariation(
        "Magnetostatic", setup_name, [], [], list(expressions),
    )

    # Magnetostatic: no time sweep, just single values
    # The sweep values might be empty or single point
    sweep_vals = np.array(list(result.GetSweepValues()), dtype=float)
    data = {}
    for expr in expressions:
        real_vals = np.array(list(result.GetRealDataValues(expr)), dtype=float)
        data[expr] = real_vals

    desktop.CloseProject(project_name)
    return sweep_vals, data


# --- CLI entry point for testing ---
if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')

    import argparse
    ap = argparse.ArgumentParser(description='COM-based AEDT data extractor')
    ap.add_argument('project', help='Project directory path')
    ap.add_argument('--design', default='5_Partial_motor_TR')
    ap.add_argument('--setup', default=None)
    ap.add_argument('--type', default='transient', choices=['transient', 'magnetostatic'])
    ap.add_argument('--expressions', nargs='*', default=None)
    ap.add_argument('--csv', default=None, help='Save to CSV')
    ap.add_argument('--fft', action='store_true', help='Run FFT on extracted data')
    ap.add_argument('--speed', type=int, default=3000, help='Rated speed for FFT')
    ap.add_argument('--pole-pairs', type=int, default=4)
    ap.add_argument('--plot', default=None, help='Save plot to PNG')
    args = ap.parse_args()

    if args.type == 'transient':
        time_s, data = extract_transient(
            args.project, args.design, args.expressions, args.setup
        )
        print(f'Extracted {len(time_s)} time points, {len(data)} channels')
        for name, vals in data.items():
            print(f'  {name}: [{vals.min():.3f}, {vals.max():.3f}]')

        if args.csv:
            with open(args.csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Time_s'] + list(data.keys()))
                for i in range(len(time_s)):
                    writer.writerow([time_s[i]] + [data[k][i] for k in data])
            print(f'Saved: {args.csv}')

        if args.fft:
            time_steps_per_cycle = 50
            freq = args.speed / 60 * args.pole_pairs
            time_step = 1 / (freq * time_steps_per_cycle)

            for name, vals in data.items():
                fft_vals = np.fft.rfft(vals)
                fft_mag = np.abs(fft_vals)
                freqs = np.fft.rfftfreq(len(vals), d=time_step)
                fi = np.argmax(fft_mag[1:]) + 1
                fr = fft_mag[fi]
                fm = fr * 2 / len(vals)
                mask = np.ones(len(fft_mag), dtype=bool)
                mask[0] = False; mask[fi] = False
                thd = np.sqrt(np.sum(fft_mag[mask]**2)) / fr * 100 if fr > 0 else 0
                print(f'{name}: fundamental={fm:.2f} V, THD={thd:.2f}%')

            ke_line = data[list(data.keys())[0]].max() * np.sqrt(3)  # rough estimate
            # Better: use fundamental
            fft_a = np.fft.rfft(data[list(data.keys())[0]])
            fft_a_mag = np.abs(fft_a)
            fi_a = np.argmax(fft_a_mag[1:]) + 1
            fund_a = fft_a_mag[fi_a] * 2 / len(data[list(data.keys())[0]])
            ke_line = fund_a / (args.speed / 1000) * np.sqrt(3)
            print(f'Ke (line) = {ke_line:.4f} V/(krpm)')
    else:
        sw, data = extract_magnetostatic(
            args.project, args.design, args.expressions, args.setup
        )
        print(f'Extracted {len(sw)} sweep points, {len(data)} channels')
