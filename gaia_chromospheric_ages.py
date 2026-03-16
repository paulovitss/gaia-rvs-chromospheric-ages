#
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import simpson
import joblib
import os


############################################ FUNCTIONS ############################################################################################

def load_spectrum(spectrum_file): # load spectrum
    lamb, flux = np.genfromtxt("spectra/"+spectrum_file, delimiter=',', unpack=True, skip_header=1, missing_values='', filling_values=np.nan)
    mask = np.isfinite(flux)
    return lamb[mask], flux[mask]
    
def check_parameter_domain(teff, feh, logg, mass): # check stellar parameters domain
    if not (4500 <= teff <= 6500):
        return False
    if not (-1.0 <= feh <= 0.4):
        return False
    if logg <= 3.0:
        return False
    if not (0.7 <= mass <= 1.4):
        return False
    return True    
    
def integral_apparent_flux(lamb, flux, line): # integrate the apparent flux on the spectrum
    delta_lambda = 1.2
    mask = np.isfinite(flux)
    lamb = lamb[mask]
    flux = flux[mask]
 
    interp_flux = interp1d(lamb, flux, kind='linear')
    lamb_fine = np.arange(line-delta_lambda, line+delta_lambda + 0.01, 0.01)
    flux_fine = interp_flux(lamb_fine)
    
    integral = simpson(flux_fine, x = lamb_fine)
    return(integral)
    
 
# Coeficients to convert the observed apparent flux to the scale of Souza dos Santos+2026
CONVERSION = {
    8498.02: (1.23013299, -0.43224595),
    8542.09: (1.29793531, -0.43653734),
    8662.14: (1.25313884, -0.4023392)}
    
def convert_flux(line, flux): # convert the flux scale
    a, b = CONVERSION[line]
    return a*flux + b
    
def interpolar_fluxos(teff, logg, feh): # calculate the absolute fluxes on the continum windows 
    ponto = np.array([[teff, logg, feh]])
    resultados = {}
    for nome, interpolador in interpolators.items():
        resultado = interpolador(ponto)
        resultados[nome] = resultado[0] if resultado.size > 0 else np.nan
    return resultados

def calcular_fluxo_absoluto(apparent_flux, teff, logg, feh, line): # convert the apparent flux to absolute flux
    fluxes = interpolar_fluxos(teff, logg, feh)
    y = np.array([fluxes["R1"], fluxes["R4"], fluxes["R6"], fluxes["R7"], fluxes["R8"]])

    coef, cov = np.polyfit(mean_lambda, y, 1, cov=True) # fitting the continuum inclination
    a, b = coef
    
    y_pred = a * mean_lambda + b
    rmse = np.sqrt(np.mean((y - y_pred)**2)) # estimating the RMSE of the continuum fit (not used right now)

    absolute_total_flux = apparent_flux*((a*line + b)) # calculating the absolute total flux
    return(round(absolute_total_flux), round(rmse))
    
def Fphot(line, teff, feh): #  photospheric corrections
    if line == 8498.02:
        env_teff = (1.069)*teff**2 + (-7.948e3)*teff + (1.695e7)
        env_feh = (1.679e5)*feh**2 + (-4.851e5)*feh + (1.124e5)
    elif line == 8542.09:
        env_teff = (1.084)*teff**2 + (-9.081e3)*teff + (2.098e7)
        env_feh = (3.466e5)*feh**2 + (-3.519e5)*feh + (-2.974e4)
    elif line == 8662.14:
        env_teff = (0.855)*teff**2 + (-6.523e3)*teff + (1.417e7)
        env_feh = (3.450e5)*feh**2 + (-3.208e5)*feh + (2.638e4)
    return env_teff + env_feh


def predict_age(flux, mass, feh, radius): # age calibration
    if flux < 5.3e5:
        return np.nan
 
    log_flux_sun = np.log10(747064.450077) # solar flux to scale the fit

    # Coefficients
    a = 9.6895      # const
    b = -3.8237     # log(M)
    c = 0.2077      # [Fe/H]
    d = 0.9767      # log(R)
    e = -1.2947     # F_TRIP
    f = -5.9711     # F_TRIP^2

    delta_logF = np.log10(flux) - log_flux_sun
    log_age = (
        a
        + b * np.log10(mass)
        + c * feh
        + d * np.log10(radius)
        + e * delta_logF
        + f * delta_logF**2)
    return log_age
        
def chromospheric_age(spectrum, teff, feh, logg, mass, radius): # chromospheric flux and age calculation
    results_lines = []
    chrom_fluxes = []
    lamb, flux = load_spectrum(spectrum)
    
    for line in LINES:
        apparent = integral_apparent_flux(lamb, flux, line)
        apparent_scaled = convert_flux(line, apparent)
        absolute_flux, _ = calcular_fluxo_absoluto(apparent_scaled, teff, logg, feh, line)
        fphot = Fphot(line, teff, feh)

        chrom_flux = absolute_flux - fphot
        chrom_fluxes.append(chrom_flux)

        results_lines.append({
            "line": line,
            "apparent": apparent,
            "scaled": apparent_scaled,
            "Ftotal": absolute_flux,
            "Fphot": fphot,
            "Fchrom": chrom_flux
        })

    chrom_flux_mean = np.nanmean(chrom_fluxes)
    age = predict_age(chrom_flux_mean, mass, feh, radius)
    return results_lines, chrom_flux_mean, age
    
def plotting(spectrum, plot, results_lines, flux_mean, age):
    if plot:
        os.makedirs("plots", exist_ok=True)
        width_plot = 25

        lamb, flux = load_spectrum(spectrum)
        mask = np.isfinite(flux)
        lamb = lamb[mask]
        flux = flux[mask]
        interp_flux = interp1d(lamb, flux, kind='linear')

        lamb_fine1 = np.arange(LINES[0]-1.2, LINES[0]+1.2+0.01, 0.01)
        lamb_fine2 = np.arange(LINES[1]-1.2, LINES[1]+1.2+0.01, 0.01)
        lamb_fine3 = np.arange(LINES[2]-1.2, LINES[2]+1.2+0.01, 0.01)

        flux_fine1 = interp_flux(lamb_fine1)
        flux_fine2 = interp_flux(lamb_fine2)
        flux_fine3 = interp_flux(lamb_fine3)
        
        age_lin = 10**age
        age_low = 10**(age - 0.14)
        age_high = 10**(age + 0.14)

        fig, axes = plt.subplots(3, 1, figsize=(7,5.5))

        axes[0].plot(lamb, flux)
        axes[0].fill_between(lamb_fine1, flux_fine1, alpha=0.4, color='orange')
        axes[0].axvline(LINES[0]-1.2, color='k', linestyle=':')
        axes[0].axvline(LINES[0]+1.2, color='k', linestyle=':')

        axes[1].plot(lamb, flux)
        axes[1].fill_between(lamb_fine2, flux_fine2, alpha=0.4, color='orange')
        axes[1].axvline(LINES[1]-1.2, color='k', linestyle=':')
        axes[1].axvline(LINES[1]+1.2, color='k', linestyle=':')
        
        axes[2].plot(lamb, flux)
        axes[2].fill_between(lamb_fine3, flux_fine3, alpha=0.4, color='orange')
        axes[2].axvline(LINES[2]-1.2, color='k', linestyle=':')
        axes[2].axvline(LINES[2]+1.2, color='k', linestyle=':')
        
        for i, ax in enumerate(axes):
            ax.set_xlim(LINES[i]-width_plot, LINES[i]+width_plot)
            ax.set_ylim(0.1, 1.05)
            ax.set_xlabel(r'$\lambda$ ($\AA$)')
            ax.set_ylabel('Normalized flux')
        for i, ax in enumerate(axes):
            Fap   = results_lines[i]["apparent"]
            Ftot  = results_lines[i]["Ftotal"]
            Fphot = results_lines[i]["Fphot"]
            Fchr  = results_lines[i]["Fchrom"]
            text = (
                f"$F_{{\\rm ap}}$ = {Fap:.3f}\n"
                f"$F_{{\\rm tot}}$ = {Ftot:.2e}\n"
                f"$F_{{\\rm phot}}$ = {Fphot:.2e}\n"
                f"$F_{{\\rm chrom}}$ = {Fchr:.2e}"
            )
            ax.text(0.03, 0.7, text,
                transform=ax.transAxes, fontsize=9, verticalalignment="top", bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))      
                
            summary = (spectrum.replace('.txt','')+'\n'
                f"$\\langle F_{{\\rm chrom}} \\rangle$ = {flux_mean:.2e}\n"
                f"$\\log\\,{{\\rm Age}}$ = {age:.2f} $\\pm$ 0.14\n"
                f"${{\\rm Age}}$ = {age_lin/1e9:.2f} "
                f"$^{{+{(age_high-age_lin)/1e9:.2f}}}_"
                f"{{-{(age_lin-age_low)/1e9:.2f}}}$ Gyr")
            fig.text(
                0.75, 0.54,
                summary,
                fontsize=10,
                va="center",
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="black"))
            
        plt.tight_layout(rect=[0,0,0.75,1])
        plt.savefig(f"plots/{spectrum.replace('.txt','')}.pdf", dpi=300, bbox_inches='tight')
        plt.close(fig)

        
############################################ RUNNING ##############################################################################################


interpolators, mean_lambda = joblib.load("interpolator_marcs_caIItriplet.pkl").values() # load the MARCS flux interpolator
LINES = [8498.02, 8542.09, 8662.14]

def main():
    parser = argparse.ArgumentParser(description="Compute chromospheric ages from Gaia spectra")
    parser.add_argument("input_file", help="Input table with stellar parameters")
    parser.add_argument("output_file", help="Output file with results")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--save", action="store_true")

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    plot = args.plot
    save = args.save

    # reading the input table:
    stars = []
    with open(input_file) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("#") or line.strip() == "":
            continue
        spectrum, teff, feh, logg, mass, radius = line.split()
        stars.append({
            "spectrum": spectrum,
            "teff": float(teff),
            "feh": float(feh),
            "logg": float(logg),
            "mass": float(mass),
            "radius": float(radius)})

    # running the calculations:
    results = []
    for star in stars:
        valid = check_parameter_domain( star["teff"], star["feh"], star["logg"], star["mass"])
        if not valid:
            print(f"{star['spectrum']} \n outside calibration domain\n")
            results.append({
            "spectrum": star["spectrum"],
            "teff": star["teff"],
            "feh": star["feh"],
            "logg": star["logg"],
            "mass": star["mass"],
            "radius": star["radius"],
            "apparent": [np.nan, np.nan, np.nan],
            "scaled": [np.nan, np.nan, np.nan],
            "ftotal": [np.nan, np.nan, np.nan],
            "fphot": [np.nan, np.nan, np.nan],
            "flux_mean": np.nan,
            "age": np.nan,
            "flag": "OUT_OF_DOMAIN"
            })
            continue
    
        results_lines, flux_mean, age = chromospheric_age(star["spectrum"], star["teff"], star["feh"], star["logg"], star["mass"], star["radius"])
 
        apparent_vals = [r["apparent"] for r in results_lines]
        scaled_vals = [r["scaled"] for r in results_lines]
        ftotal_vals = [r["Ftotal"] for r in results_lines]
        fphot_vals = [r["Fphot"] for r in results_lines]
        fchrom_vals = [r["Fchrom"] for r in results_lines]
        
        plotting(star['spectrum'], plot, results_lines, flux_mean, age)
        
        if flux_mean < 5.3e5:
            flag = 'OUT_OF_DOMAIN'
        elif (flux_mean < 6e5)*(flux_mean > 5.3e5):
            flag = 'LOW_ACTIVITY_REGIME'
        else:
            flag = 'OK'

        results.append({
            "spectrum": star["spectrum"],
            "teff": star["teff"],
            "feh": star["feh"],
            "logg": star["logg"],
            "mass": star["mass"],
            "radius": star["radius"],
            "apparent": apparent_vals,
            "scaled": scaled_vals,
            "ftotal": ftotal_vals,
            "fphot": fphot_vals,
            "fchrom": fchrom_vals,
            "flux_mean": flux_mean,
            "age": age,
            "flag": flag
            })

        print(star["spectrum"])
        print(f"Chromospheric flux: {flux_mean:.3e}")
        if flux_mean < 5.3e5:
            print("Warning: chromospheric flux below the age calibration limit.\n")
        else: 
            print(f"log Age: {age:.3f} ± 0.14 dex\n")

    # writing output file:
    with open(output_file, "w") as f:
        f.write("spectrum teff feh logg mass rad "
        "Fap_8498 Fap_8542 Fap_8662 "
        "Fsc_8498 Fsc_8542 Fsc_8662 "
        "Ftot_8498 Ftot_8542 Ftot_8662 "
        "Fphot_8498 Fphot_8542 Fphot_8662 "
        "Fchrom_8498 Fchrom_8542 Fchrom_8662 "
        "Fchrom_mean log_age flag\n")
        
        for r in results:
            apparent_vals = r["apparent"]
            scaled_vals = r["scaled"]
            ftotal_vals = r["ftotal"]
            fphot_vals = r["fphot"]
            fchrom_vals = r["fchrom"]

            f.write(
                f"{r['spectrum']} "
                f"{r['teff']:.0f} {r['feh']:.2f} {r['logg']:.2f} {r['mass']:.2f} {r['radius']:.2f} "
                f"{apparent_vals[0]:.3f} {apparent_vals[1]:.3f} {apparent_vals[2]:.3f} "
                f"{scaled_vals[0]:.3f} {scaled_vals[1]:.3f} {scaled_vals[2]:.3f} "
                f"{ftotal_vals[0]:.3e} {ftotal_vals[1]:.3e} {ftotal_vals[2]:.3e} "
                f"{fphot_vals[0]:.3e} {fphot_vals[1]:.3e} {fphot_vals[2]:.3e} "
                f"{fchrom_vals[0]:.3e} {fchrom_vals[1]:.3e} {fchrom_vals[2]:.3e} "
                f"{r['flux_mean']:.3e} {r['age']:.3f} {r['flag']}\n"
                )

if __name__ == "__main__":
    main()


