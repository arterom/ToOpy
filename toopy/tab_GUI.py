import tkinter
import tkinter.messagebox
import customtkinter
import subprocess
from PIL import Image
import os
import time
from tkinter.messagebox import askokcancel, showinfo, WARNING
customtkinter.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # configure window
        self.title("GUI for Reference Alerts.py")
        self.geometry(f"{1000}x{1000}")

        # configure grid layout (4x4)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure((0, 1, 2), weight=1)

        self.textbox = customtkinter.CTkTextbox(self, width=1000, height=200)
        self.textbox.grid(row=0, column=0, padx=(20, 0), pady=(20, 0), sticky="nsew")

        self.tabview = customtkinter.CTkTabview(self, width=1000)
        self.tabview.grid(row=1, column=0, padx=(20, 20), pady=(20, 0), sticky="nsew")
        #self.tabview = customtkinter.CTkLabel(master=self.tabview, text="ghjgjhgj Group:")
        self.tabview.add("GW Alert")
        self.tabview.add("Neutrino Alert")
        self.tabview.add("Swift-BAT")
        self.tabview.add("Fermi-GBM")



        #self.string_input_button_swift = customtkinter.CTkButton(self.tabview.tab("Swift-BAT"), command=self.random, text="Execute for Swift-BAT")
        #self.string_input_button_swift.grid(row=1, column=0, padx=20, pady=(10, 10))





        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./GUI_images")
        self.large_test_image = customtkinter.CTkImage(Image.open(os.path.join(image_path, "ToOpy.png")), size=(700, 300))
        self.home_frame_large_image_label = customtkinter.CTkLabel(self, text="", image=self.large_test_image)
        self.home_frame_large_image_label.grid(row=2, column=0, padx=(20, 0), pady=(20, 0), sticky="nsew")


        self.string_input_button = customtkinter.CTkButton(self.tabview.tab("GW Alert"), command=self.input_gw_event, text="Execute for GW")
        self.string_input_button.grid(row=1, column=0, padx=20, pady=(10, 10))
        self.button_2 = customtkinter.CTkButton(self.tabview.tab("GW Alert"), command=self.initalize_Reference_GW, text="Run GW A")
        self.button_2.grid(row=2, column=0, padx=20, pady=(10, 10))



        #self.label_tab_2 = customtkinter.CTkLabel(self.tabview.tab("Neutrino Alert"), text="IceCube Track")
        self.string_input_button_Track = customtkinter.CTkButton(self.tabview.tab("Neutrino Alert"), text="Execute for Track",
                                                           command=self.input_neutrino_event_Track)
        self.string_input_button_Track.grid(row=0, column=0, padx=20, pady=(10, 10))
        #self.label_tab_2.grid(row=0, column=0, padx=20, pady=20)


        #self.label_tab_2 = customtkinter.CTkLabel(self.tabview.tab("Neutrino Alert"), text="IceCube Cascade")
        self.string_input_button_Cascade = customtkinter.CTkButton(self.tabview.tab("Neutrino Alert"), text="Execute for IC170922A_HE",
                                                           command=self.input_neutrino_event_IC170922A_HE)
        self.string_input_button_Cascade.grid(row=1, column=0, padx=20, pady=(10, 10))
        #self.label_tab_2.grid(row=2, column=0, padx=20, pady=20)

        self.button_2 = customtkinter.CTkButton(self.tabview.tab("Neutrino Alert"), command=self.initalize_Reference_IC170922A_ROI, text="Run Track VarInd")
        self.button_2.grid(row=2, column=0, padx=20, pady=(10, 10))

        self.textbox.insert("0.0","Experimental GUI in order to facilitate the use of ToOpy (for now works with reference alerts that are processed through .sh scripts)." "\n\n" + 
        	"How to Use:\n" + "1). Select alert streams from the tabview below and run pre-defined scripts by pushing run buttons\n" + 
        	"2). Process reference alerts by pushing the corresponding buttons \n" + "3). Check updated output folders shown in the figure at the bottom (note that ToOpy might be still running (especially for Fermi-LAT analysis)")


    def initalize_Reference_GW(self):
        print('x')
        subprocess.check_call('chmod u+r+x initalize_Reference_GW_A.sh', shell=True)
        subprocess.check_call('./initalize_Reference_GW_A.sh', shell=True)
        print('at the end')
        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./GUI_images")
        self.large_test_image = customtkinter.CTkImage(Image.open(os.path.join(image_path, "ToOpy.png")), size=(1000, 300))
        self.home_frame_large_image_label = customtkinter.CTkLabel(self, text="", image=self.large_test_image)
        self.home_frame_large_image_label.grid(row=2, column=0, padx=(20, 0), pady=(20, 0), sticky="nsew")
    

    def initalize_Reference_IC170922A_ROI(self):
        print('x')
        subprocess.check_call('chmod u+r+x initalize_Reference_IC170922A_ROI.sh', shell=True)
        subprocess.check_call('./initalize_Reference_IC170922A_ROI.sh', shell=True)
        print('at the end')
        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./GUI_images")
        self.large_test_image = customtkinter.CTkImage(Image.open(os.path.join(image_path, "track.png")), size=(1000, 300))
        self.home_frame_large_image_label = customtkinter.CTkLabel(self, text="", image=self.large_test_image)
        self.home_frame_large_image_label.grid(row=2, column=0, padx=(20, 0), pady=(20, 0), sticky="nsew")

        file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./")
        for i in range(0,2):
        	print(i)



    def input_gw_event(self):
        dialog = customtkinter.CTkInputDialog(text="Add conditions in the following format: ...", title="Confirmation")
        print("This is fed into .sh Scripts:", dialog.get_input())

    #def input_neutrino_event_IC170922A_HE(self):
    #    dialog = customtkinter.CTkInputDialog(text="Carefull this will take time to run", title="Progress")
    #    print("This is fed into .sh Scripts:", dialog.get_input())
    #    subprocess.check_call('chmod u+r+x initalize_Reference_IC170922A_HE.sh', shell=True)
    #    subprocess.check_call('./initalize_Reference_IC170922A_HE.sh', shell=True)

    def input_neutrino_event_IC170922A_HE(self):
        print('x')
        answer = askokcancel(
            title='Confirmation',
            message='Running Fermi-LAT Analysis will take some time.',
            icon=WARNING)

        if answer:
            print('yes to running')
            subprocess.check_call('chmod u+r+x initalize_Reference_IC170922A_HE.sh', shell=True)
            subprocess.check_call('./initalize_Reference_IC170922A_HE.sh', shell=True)
            showinfo(
                title='Deletion Status',
                message='ToOpy is processing the alert please check command line and output folders')

    def input_neutrino_event_Track(self):
        dialog = customtkinter.CTkInputDialog(text="VarInd or HE?", title="Progress")
        print("This is fed into .sh Scripts:", dialog.get_input())
        if (print(dialog.get_input()) == "VarInd"):
        	print('x')



        argument_list = dialog.get_input()
        print(argument_list)
        #argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skymap_fits_url), 'VarInd', fermitools_refdata_path, 'no']
        separator = " "
        subprocess.check_call("./method_scripts/IceCube_TRACK.sh %s" % separator.join(argument_list), shell=True)

if __name__ == "__main__":
    app = App()
    app.mainloop()





