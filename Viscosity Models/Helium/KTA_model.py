import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from helium_visc_data import viscData

#from https://www.sciencedirect.com/science/article/pii/S0149197024004670#bib28:
def KTA(T):
    eta = 3.674*10**(-7) * T**(0.7)    #[Pa*s]
    return eta

#from https://www.sciencedirect.com/science/article/pii/S0149197024004670#bib28:
def KTA_mod(T):
    eta = 3.817*10**(-7) * T**(0.6938)   #[Pa*s]
    return eta

### DETTE ER DEN NYE VISKOSITETSMODELLEN FOR HELIUM ###
def KTA_tweak(T, P):
    P_crit = 0.22832   #[MPa]   (Source: NIST)
    #eta = 1e-7*(3.817 * T**(0.6938) + (P**((2 - T/300)**(5.05)) / (T*P_crit)) - (np.exp(-(T-325)**2 / 1000))*(T**((2-300/T)**(2) - 1)) + (np.exp(-(T-325)**2 / 1000))*(P/25)**(2.7)) * 10**6   #[µPa*s]
    A = (2 - T/300)**(5.05)
    B = (2-300/T)**(2) - 1
    eta = 1e-7 * (3.817 * T**(0.6938) + P**A/(T*P_crit) + (np.exp(-(T-325)**2 / 1000))*((P/25)**(2.7) - T**B ))   #[Pa*s]
    return eta

def getModelViscosity(T, P, viscosityModel):
    if viscosityModel == "KTA":
        eta = KTA(T) * 10**6    #[µPa*s]
    elif viscosityModel == "KTA_mod":
        eta = KTA_mod(T) * 10**6    #[µPa*s]
    elif viscosityModel == "KTA_tweak":
        eta = KTA_tweak(T, P) * 10**6    #[µPa*s]
    else:
        raise ValueError("Invalid viscosity model")
    
    return eta

def modelPerformance(T, viscosityModels):
    #dataList object: [Temperature, Pressure, Viscosity, Author, ModelViscosity, ARD, markerSymbol (for plot)]
    for viscosityModel in viscosityModels:
        dataList = ARDViscIsotherm(T, viscosityModel)
        if not dataList:
            print(f"No data found at {T}K")
            break
        averageARD = np.mean([abs(row[5]) for row in dataList])
        maxARD = max([abs(row[5]) for row in dataList])
        minPressure = min([row[1] for row in dataList])
        maxPressure = max([row[1] for row in dataList])
        print(f"Model: {viscosityModel}")
        print(f"Temperature: {T:.2f}K || Pressure range: {minPressure:.2f}MPa - {maxPressure:.2f}MPa")
        print(f"Average ARD: {averageARD:.2f}%")
        print(f"Max ARD: {maxARD:.2f}% \n")


def ARDViscIsotherm(targetTemp, viscModel):
    data_dict = viscData()
    dataList = []
    for key, data in data_dict.items():
        for i in range(len(data[0])):
            margin = abs(data[0][i] - targetTemp)
            if margin <= 0.5:
                #data_row = (data[0].iloc[i], data[1].iloc[i], data[2].iloc[i], key)
                modelViscosity = getModelViscosity(data[0].iloc[i], data[1].iloc[i], viscModel)
                ARD = 100 * (modelViscosity - data[2].iloc[i]) / data[2].iloc[i]
                data_row = [data[0].iloc[i], data[1].iloc[i], data[2].iloc[i], key, modelViscosity, ARD, data[3]]
                #data_row object: [Temperature, Pressure, Viscosity, Author, ModelViscosity, ARD, markerSymbol (for plot)]
                dataList.append(data_row)

    return dataList



def ARDPlot(Temperatures, viscModels):
    fig, axs = plt.subplots(2, 2, figsize=(14, 8))
    axs = axs.flatten()

    # Define colors for each viscosity model
    colors = {
        "KTA": "black",
        "KTA_mod": "red",
        "KTA_tweak": "blue"
    }

    # Dictionary to keep track of labels added to the legend
    legend_labels = {}

    for i, T in enumerate(Temperatures):
        for viscModel in viscModels:
            dataList = ARDViscIsotherm(T, viscModel)
            for row in dataList:
                pressures = row[1]
                ARDs = row[5]
                markerSymbol = row[6]
                facecolor = colors[viscModel] if markerSymbol in ['p', 'x', '2', '+'] else 'none'
                label = f"{viscModel} ({row[3]})"
                if label not in legend_labels:
                    axs[i].scatter(pressures, ARDs, label=label, marker=markerSymbol, edgecolors=colors[viscModel], facecolors=facecolor)
                    legend_labels[label] = axs[i].scatter([], [], label=label, marker=markerSymbol, edgecolors=colors[viscModel], facecolors=facecolor)
                else:
                    axs[i].scatter(pressures, ARDs, marker=markerSymbol, edgecolors=colors[viscModel], facecolors=facecolor)
            axs[i].axhline(y=0, color='black', linestyle='-', linewidth=0.7)
        axs[i].set_ylim(-1.5, 1.5)  # Set y-axis limits
        if i >= 2:  # Only set x-label for the bottom plots
            axs[i].set_xlabel("Pressure [MPa]")
        axs[i].set_title(f"{T}K")
        axs[i].grid()
    
    fig.text(0, 0.5, r'$100 \left(\eta_{\mathrm{model}} - \eta_{\mathrm{experimental}}\right)/\eta_{\mathrm{experimental}}$', 
             ha='center', va='center', rotation='vertical', fontsize=15)    
    fig.suptitle("ARD Plots", fontsize=16, y=0.93)

    # Create a single legend from the plot elements
    handles, labels = [], []
    for handle, label in legend_labels.items():
        handles.append(label)
        labels.append(handle)
    fig.legend(handles, labels, loc='lower center', ncol=5, bbox_to_anchor=(0.5, 0.0), frameon=False)

    fig.tight_layout(rect=[0, 0.1, 1, 0.95])
    plt.savefig("helium_grouped_plot.png", dpi=300, bbox_inches='tight')
    plt.clf()



ARDPlot([273.15, 293.15, 323.15, 373.15], ["KTA", "KTA_mod", "KTA_tweak"])

modelPerformance(293.15, ["KTA", "KTA_mod", "KTA_tweak"])
