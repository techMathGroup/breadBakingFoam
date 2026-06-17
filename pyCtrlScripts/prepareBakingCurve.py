# -- python script to create baking curve according to parameters
# -- INTERACTIVE VERSION with sliders

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import os

# rozume pars = [260 - 300], [5 - 15], [30 - 70]

def createBakingCurve(pars, bakingLength=40, nPoints=200):
    """Calculate baking curve from parameters"""
    initTemp, lenInit, decayRate = pars
    
    time = np.linspace(0, bakingLength, nPoints)
    bakingCurve = initTemp + 0*time

    for i, t in enumerate(time):
        if t > lenInit:
            bakingCurve[i] = initTemp - decayRate * np.log10((t - lenInit) + 1)        

    return time, bakingCurve


def interactive_baking_curve():
    """Interactive plot with parameter sliders"""
    
    # Initial parameters
    init_initTemp = 280
    init_lenInit = 10
    init_decayRate = 70
    bakingLength = 40
    
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(12, 8))
    plt.subplots_adjust(left=0.1, bottom=0.35, right=0.9, top=0.95)
    
    # Initial plot
    time, bakingCurve = createBakingCurve([init_initTemp, init_lenInit, init_decayRate], bakingLength, nPoints=200)
    line, = ax.plot(time, bakingCurve, 'r-', linewidth=2.5, label='Baking Temperature')
    
    # Add reference lines
    ax.axhline(y=100, color='blue', linestyle='--', alpha=0.3, label='Water boiling (100°C)')
    ax.axhline(y=180, color='green', linestyle='--', alpha=0.3, label='Typical final temp (180°C)')
    
    ax.set_xlabel('Baking Time (min)', fontsize=12)
    ax.set_ylabel('Temperature (°C)', fontsize=12)
    ax.set_title('Interactive Baking Temperature Profile', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    ax.set_xlim([0, bakingLength])
    ax.set_ylim([130, 300])
    
    # Create slider axes
    ax_initTemp = plt.axes([0.15, 0.20, 0.7, 0.03])
    ax_lenInit = plt.axes([0.15, 0.15, 0.7, 0.03])
    ax_decayRate = plt.axes([0.15, 0.10, 0.7, 0.03])
    
    # Create sliders
    slider_initTemp = Slider(
        ax=ax_initTemp,
        label='Initial Temp (°C)',
        valmin=200,
        valmax=320,
        valinit=init_initTemp,
        valstep=5,
        color='orangered'
    )
    
    slider_lenInit = Slider(
        ax=ax_lenInit,
        label='Hold Time (min)',
        valmin=0,
        valmax=20,
        valinit=init_lenInit,
        valstep=1,
        color='steelblue'
    )
    
    slider_decayRate = Slider(
        ax=ax_decayRate,
        label='Decay Rate',
        valmin=20,
        valmax=120,
        valinit=init_decayRate,
        valstep=5,
        color='forestgreen'
    )
    
    # Update function
    def update(val):
        initTemp = slider_initTemp.val
        lenInit = slider_lenInit.val
        decayRate = slider_decayRate.val
        
        time, bakingCurve = createBakingCurve([initTemp, lenInit, decayRate], bakingLength)
        line.set_ydata(bakingCurve)
        
        # Update title with current values
        ax.set_title(f'Baking Profile: {initTemp:.0f}°C hold for {lenInit:.0f} min, then decay with rate {decayRate:.0f}',
                    fontsize=11, fontweight='bold')
        
        fig.canvas.draw_idle()
    
    # Connect sliders to update function
    slider_initTemp.on_changed(update)
    slider_lenInit.on_changed(update)
    slider_decayRate.on_changed(update)
    
    # Add save button
    ax_save = plt.axes([0.4, 0.03, 0.2, 0.04])
    btn_save = Button(ax_save, 'Save Curve', color='lightgoldenrodyellow', hovercolor='gold')
    
    def save_curve(event):
        initTemp = slider_initTemp.val
        lenInit = slider_lenInit.val
        decayRate = slider_decayRate.val
        
        # Create directory if needed
        os.makedirs('bakingCurves', exist_ok=True)
        
        # Save figure
        fig.savefig('bakingCurves/interactive_curve.png', dpi=150, bbox_inches='tight')
        
        # Save data to CSV
        time, bakingCurve = createBakingCurve([initTemp, lenInit, decayRate], bakingLength)
        np.savetxt('bakingCurves/curve_data.csv', 
                   np.column_stack([time, bakingCurve]),
                   delimiter=',',
                   header='time_min,temperature_C',
                   comments='')
        
        print(f"\nSaved:")
        print(f"  - Image: bakingCurves/interactive_curve.png")
        print(f"  - Data:  bakingCurves/curve_data.csv")
        print(f"  - Parameters: initTemp={initTemp:.0f}, lenInit={lenInit:.0f}, decayRate={decayRate:.0f}\n")
    
    btn_save.on_clicked(save_curve)
    
    # Add reset button
    ax_reset = plt.axes([0.7, 0.03, 0.15, 0.04])
    btn_reset = Button(ax_reset, 'Reset', color='lightcoral', hovercolor='salmon')
    
    def reset(event):
        slider_initTemp.reset()
        slider_lenInit.reset()
        slider_decayRate.reset()
    
    btn_reset.on_clicked(reset)
    
    # Add instruction text
    fig.text(0.5, 0.25, 'Adjust sliders to modify the baking temperature profile', 
             ha='center', fontsize=10, style='italic', color='gray')
    
    plt.show()


if __name__ == "__main__":
    print("\n" + "="*60)
    print("INTERACTIVE BAKING CURVE GENERATOR")
    print("="*60)
    print("\nAdjust the sliders to modify the temperature profile:")
    print("  • Initial Temp: Starting oven temperature (°C)")
    print("  • Hold Time: Duration at initial temperature (min)")
    print("  • Decay Rate: Speed of temperature decrease")
    print("\nClick 'Save Curve' to export the current profile.")
    print("="*60 + "\n")
    
    interactive_baking_curve()