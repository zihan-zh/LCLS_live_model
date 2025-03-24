import matplotlib.pyplot as plt
import numpy as np
import torch
from scipy.interpolate import interp1d


def find_element_index_by_name(module_list, name):
    name_lower = name.lower()
    for idx, element in enumerate(module_list):
        if hasattr(element, 'name') and element.name.lower() == name_lower:
            print(f"Element '{element.name}' found at index: {idx}")
            return idx
    return -1  # If not found

def get_cheetah_beam_values(segment, beam):
    """Get twiss parameter evolution along the segment."""
    longitudinal_beams = [beam]
    s_positions = [torch.tensor([0.0])]
    ele_name = []
    for element in segment.elements:
        # if element.length == 0:
            # continue

        outgoing = element.track(longitudinal_beams[-1])
        ele_name.append(element.name)
        longitudinal_beams.append(outgoing)
        s_positions.append(s_positions[-1] + element.length)

    results = {
        's': s_positions,
        'beta_x': [beam.beta_x for beam in longitudinal_beams],
        'beta_y': [beam.beta_y for beam in longitudinal_beams],
        'alpha_x': [beam.alpha_x for beam in longitudinal_beams],
        'alpha_y': [beam.alpha_y for beam in longitudinal_beams],
        'sigma_t': [beam.sigma_tau for beam in longitudinal_beams],
        'energy': [beam.energy for beam in longitudinal_beams],
        'mu_x': [beam.mu_x for beam in longitudinal_beams],
        'mu_y': [beam.mu_y for beam in longitudinal_beams],
        'sigma_x': [beam.sigma_x for beam in longitudinal_beams],
        'sigma_y': [beam.sigma_y for beam in longitudinal_beams],
        'emit_x': [beam.emittance_x for beam in longitudinal_beams],
        'emit_y': [beam.emittance_y for beam in longitudinal_beams],
        'ele_name': ele_name
    }
    
    return results

def plot_twiss_parameters(bmad_output, cheetah_output, track_start_element_name, xlims=None, ylims=None):
    """Plots Twiss parameters from cheetah_output and bmad_output.
    
    Parameters:
    bmad_output (dict): Dictionary containing Bmad tracking output.
    cheetah_output (dict): Dictionary containing Cheetah tracking output.
    track_start_element_name (str): Name of the element where tracking starts.
    xlims (list): Optional x-axis limits for the plots.
    ylims (list): Optional y-axis limits for the plots.
    """
    
    # Ensure track_start_element_name is stripped of possible leading/trailing spaces and is uppercase
    track_start_element_name = track_start_element_name.strip().upper()
    
    # Find Bmad start index
    name_array = np.array([name.strip().upper() for name in bmad_output['ele.name']])
    bmad_start_index = np.where(name_array == track_start_element_name)[0]
    
    if len(bmad_start_index) == 0:
        raise ValueError(f"Element name '{track_start_element_name}' not found in Bmad output.")
    
    bmad_start_s = bmad_output['ele.s'][bmad_start_index[0]]
    
    s = cheetah_output['s']
    cheetah_beta_x = cheetah_output['beta_x']
    cheetah_beta_y = cheetah_output['beta_y']
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))  
    
    # Plot Cheetah output
    ax1.set_title(f'Tracking from {track_start_element_name} using Cheetah')
    ax1.plot(s + bmad_start_s, cheetah_beta_x, label=r'$\beta_x$')
    ax1.plot(s + bmad_start_s, cheetah_beta_y, label=r'$\beta_y$')
    ax1.set_ylabel('Twiss Beta X/Y [m]')
    
    if xlims:
        ax1.set_xlim(xlims)
    if ylims:
        ax1.set_ylim(ylims)
    
    ax1.grid()
    ax1.legend()
    
    # Plot Bmad output
    ax2.set_title(f'Tracking from {track_start_element_name} using Bmad')
    ax2.set_xlabel('s [m]')
    ax2.set_ylabel('Twiss Beta X/Y [m]')
    ax2.plot(bmad_output['ele.s'], bmad_output['ele.a.beta'], label=r'$\beta_x$')
    ax2.plot(bmad_output['ele.s'], bmad_output['ele.b.beta'], label=r'$\beta_y$')
    
    if xlims:
        ax2.set_xlim(xlims)
    if ylims:
        ax2.set_ylim(ylims)
    
    ax2.grid()
    ax2.legend()
    
    plt.show()


def plot_twiss_comparison(bmad_output, cheetah_output, track_start_element_name, track_end_element_name, locations=None, xlims=None, ylims=None, energy_lims=None):
    """Plots Twiss parameters from cheetah_output and bmad_output with optional location markers and axis limits.
    
    Parameters:
    bmad_output (dict): Dictionary containing Bmad tracking output.
    cheetah_output (dict): Dictionary containing Cheetah tracking output.
    track_start_element_name (str): Name of the element where tracking starts.
    locations (list, optional): List of tuples with element name and location for marking (element, location).
    xlims (list, optional): Optional x-axis limits for the plots.
    ylims (list, optional): Optional y-axis limits for the Twiss parameter plots.
    energy_lims (list, optional): Optional y-axis limits for the Beam energy plot.
    """
    
    # Ensure track_start_element_name is stripped of possible leading/trailing spaces and is uppercase
    track_start_element_name = track_start_element_name.strip().upper()
    
    # Find Bmad start index
    name_array = np.array([name.strip().upper() for name in bmad_output['ele.name']])
    bmad_start_index = np.where(name_array == track_start_element_name.upper())[0]
    bmad_end_index = np.where(name_array == track_end_element_name.upper())[0]
    
    if len(bmad_start_index) == 0:
        raise ValueError(f"Element name '{track_start_element_name}' not found in Bmad output.")
    
    bmad_start_s = bmad_output['ele.s'][bmad_start_index[0]]
    
    s = cheetah_output['s']
    cheetah_beta_x = cheetah_output['beta_x']
    cheetah_beta_y = cheetah_output['beta_y']
    cheetah_energy = cheetah_output['energy']
    
    # Interpolation and deviation calculation
    s_np = np.array([tensor.cpu().item() if tensor.numel() == 1 else tensor.cpu().numpy() for tensor in s])

    # Define a common set of points based on x-axis limits
    common_s = np.linspace(bmad_output['ele.s'][bmad_start_index[0]], bmad_output['ele.s'][bmad_end_index[0]], num=3000)

    # Interpolate the curves
    cheetah_interp_x = interp1d(s_np + bmad_start_s, np.array(cheetah_beta_x).flatten(), kind='linear', fill_value='extrapolate')
    bmad_interp_x = interp1d(bmad_output['ele.s'], bmad_output['ele.a.beta'], kind='linear', fill_value='extrapolate')

    cheetah_interp_y = interp1d(s_np + bmad_start_s, np.array(cheetah_beta_y).flatten(), kind='linear', fill_value='extrapolate')
    bmad_interp_y = interp1d(bmad_output['ele.s'], bmad_output['ele.b.beta'], kind='linear', fill_value='extrapolate')

    # Calculate the difference at common points
    difference_x = cheetah_interp_x(common_s) - bmad_interp_x(common_s)
    difference_y = cheetah_interp_y(common_s) - bmad_interp_y(common_s)

    # Plotting
    fig, axs = plt.subplots(5, 1, figsize=(10, 20))  # Include additional subplots for differences
    ax1, ax2, ax3, ax4, ax5 = axs
    
    # Plot Beta X
    ax1.set_title(f'Tracking from {track_start_element_name}')
    ax1.plot(bmad_output['ele.s'], bmad_output['ele.a.beta'], linestyle='--', label=r'Bmad $\beta_x$')
    ax1.plot(s + bmad_start_s, cheetah_beta_x, linestyle='--', label=r'Cheetah $\beta_x$')
    ax1.set_xlabel('s [m]')
    ax1.set_ylabel('Twiss Beta X [m]')
    
    if xlims:
        ax1.set_xlim(xlims)
    if ylims:
        ax1.set_ylim(ylims)
    
    ax1.grid()
    ax1.legend()
    
    # Plot Beta Y
    ax2.set_xlabel('s [m]')
    ax2.set_ylabel('Twiss Beta Y [m]')
    ax2.plot(bmad_output['ele.s'], bmad_output['ele.b.beta'], linestyle='--', label=r'Bmad $\beta_y$')
    ax2.plot(s + bmad_start_s, cheetah_beta_y, linestyle='--', label=r'Cheetah $\beta_y$')
    
    if xlims:
        ax2.set_xlim(xlims)
    if ylims:
        ax2.set_ylim(ylims)
    
    ax2.grid()
    ax2.legend()

    # Plot Beam Energy
    ax3.set_xlabel('s [m]')
    ax3.set_ylabel('Beam Energy [eV]')
    ax3.plot(s + bmad_start_s, cheetah_energy, label=r'Cheetah Energy')
    ax3.plot(bmad_output['ele.s'], bmad_output['ele.p0c'], linestyle='--', label=r'Bmad Energy')
    
    if xlims:
        ax3.set_xlim(xlims)
    if energy_lims:
        ax3.set_ylim(energy_lims)
    
    ax3.grid()
    ax3.legend()
    
    # Plot Difference in Beta X
    ax4.set_title('Difference in Beta X')
    ax4.plot(common_s, difference_x, label=r'$\Delta \beta_x$')
    ax4.set_xlabel('s [m]')
    ax4.set_ylabel('Difference in Beta X [m]')
    
    if xlims:
        ax4.set_xlim(xlims)
    
    ax4.grid()
    ax4.legend()
    
    # Plot Difference in Beta Y
    ax5.set_title('Difference in Beta Y')
    ax5.plot(common_s, difference_y, label=r'$\Delta \beta_y$')
    ax5.set_xlabel('s [m]')
    ax5.set_ylabel('Difference in Beta Y [m]')
    
    if xlims:
        ax5.set_xlim(xlims)
    
    ax5.grid()
    ax5.legend()
    
    # Mark locations on Beta X, Y, and Difference plots, if provided
    if locations:
        for element, location in locations:
            if location is not None:
                ax1.axvline(x=location, color='red', linestyle='--', linewidth=1)
                ax2.axvline(x=location, color='red', linestyle='--', linewidth=1)
                ax3.axvline(x=location, color='red', linestyle='--', linewidth=1)
                ax4.axvline(x=location, color='red', linestyle='--', linewidth=1)
                ax5.axvline(x=location, color='red', linestyle='--', linewidth=1)
                ax2.text(location, ax2.get_ylim()[1] * 0.9, element, rotation=90, verticalalignment='top', color='red')
                ax5.text(location, ax5.get_ylim()[1] * 0.9, element, rotation=90, verticalalignment='top', color='red')
    
    plt.tight_layout()
    plt.show()