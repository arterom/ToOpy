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
        self.geometry(f"{900}x{1000}")

        # configure grid layout (4x4)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure((0, 1, 2), weight=1)

        self.textbox = customtkinter.CTkTextbox(self, width=800, height=200)
        self.textbox.grid(row=0, column=0, padx=(20, 0), pady=(20, 0), sticky="nsew")

        self.tabview = customtkinter.CTkTabview(self, width=800)
        self.tabview.grid(row=1, column=0, padx=(20, 20), pady=(20, 0), sticky="nsew")
        #self.tabview = customtkinter.CTkLabel(master=self.tabview, text="ghjgjhgj Group:")
        self.tabview.add("GW Alert")
        self.tabview.add("Neutrino Alert")
        self.tabview.add("Swift-BAT Alert")
        self.tabview.add("Fermi-GBM Alert")



        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./documentation_images")
        self.large_test_image = customtkinter.CTkImage(Image.open(os.path.join(image_path, "ToOpy.png")), size=(700, 300))
        self.home_frame_large_image_label = customtkinter.CTkLabel(self, text="", image=self.large_test_image)
        self.home_frame_large_image_label.grid(row=2, column=0, padx=(20, 0), pady=(20, 0), sticky="nsew")


        self.Button_GW_tilling = customtkinter.CTkButton(self.tabview.tab("GW Alert"), command=self.input_gw_event_GW_tilling, text="Run tilling-approach")
        self.Button_GW_tilling.grid(row=0, column=0, padx=20, pady=(10, 10))

        self.Button_GW_targeted = customtkinter.CTkButton(self.tabview.tab("GW Alert"), command=self.input_gw_event_GW_targeted, text="Run targeted-approach")
        self.Button_GW_targeted.grid(row=1, column=0, padx=20, pady=(10, 10))

        self.Button_GW170817_STMOC = customtkinter.CTkButton(self.tabview.tab("GW Alert"), command=self.input_gw_event_GW170817_STMOC, text="Run STMOC-approach (for GW170817)")
        self.Button_GW170817_STMOC.grid(row=2, column=0, padx=20, pady=(10, 10))



        self.Button_IC170922A_HE = customtkinter.CTkButton(self.tabview.tab("Neutrino Alert"), text="Run for IC170922A_HE",command=self.input_neutrino_event_IC170922A_HE)
        self.Button_IC170922A_HE.grid(row=0, column=0, padx=20, pady=(10, 10))

        self.Button_IC170922A_ROI = customtkinter.CTkButton(self.tabview.tab("Neutrino Alert"), text="Run for IC170922A_ROI",command=self.input_neutrino_event_IC170922A_ROI)
        self.Button_IC170922A_ROI.grid(row=1, column=0, padx=20, pady=(10, 10))


        self.Button_IC220303A_HE = customtkinter.CTkButton(self.tabview.tab("Neutrino Alert"), text="Run for IC220303A_HE (with LC)",command=self.input_neutrino_event_IC220303A_HE)
        self.Button_IC220303A_HE.grid(row=2, column=0, padx=20, pady=(10, 10))

        self.Button_IC220303A_ROI = customtkinter.CTkButton(self.tabview.tab("Neutrino Alert"), text="Run for IC220303A_ROI",command=self.input_neutrino_event_IC220303A_ROI)
        self.Button_IC220303A_ROI.grid(row=3, column=0, padx=20, pady=(10, 10))





        self.Button_GBM_targeted = customtkinter.CTkButton(self.tabview.tab("Fermi-GBM Alert"), command=self.input_event_GBM_targeted, text="Run targeted-approach")
        self.Button_GBM_targeted.grid(row=0, column=0, padx=20, pady=(10, 10))

        self.Button_GBM_STMOC = customtkinter.CTkButton(self.tabview.tab("Fermi-GBM Alert"), command=self.input_event_GBM_STMOC, text="Run STMOC-approach")
        self.Button_GBM_STMOC.grid(row=1, column=0, padx=20, pady=(10, 10))

        self.Button_BAT_targeted = customtkinter.CTkButton(self.tabview.tab("Swift-BAT Alert"), command=self.input_event_BAT_targeted, text="Run targeted-approach")
        self.Button_BAT_targeted.grid(row=0, column=0, padx=20, pady=(10, 10))

        self.Button_BAT_STMOC = customtkinter.CTkButton(self.tabview.tab("Swift-BAT Alert"), command=self.input_event_BAT_STMOC, text="Run STMOC-approach")
        self.Button_BAT_STMOC.grid(row=1, column=0, padx=20, pady=(10, 10))

        #self.button_2 = customtkinter.CTkButton(self.tabview.tab("Neutrino Alert"), command=self.initalize_Reference_IC170922A_ROI, text="Run Track VarInd")
        #self.button_2.grid(row=2, column=0, padx=20, pady=(10, 10))

        self.textbox.insert("0.0","Experimental GUI in order to facilitate the use of ToOpy (for now works with reference alerts that are processed through .sh scripts)." "\n\n" + 
        	"How to use:\n" + "1). Select alert streams from the tabview below\n" + 
        	"2). Process reference alerts from pre-defined '*.sh' scripts by pushing the corresponding buttons \n" + "3). Check command-line feedback and corresponding output folders in your home directory (note that ToOpy might be still running)")

    def input_gw_event_GW_tilling(self):
        answer = askokcancel(
            title='Confirmation',
            message='Computing tiling sequence for GW alert and checking observability from Roque de los Muchachos.',
            icon=WARNING)
        if answer:
            subprocess.check_call('chmod u+r+x initalize_Reference_GW_tilling.sh', shell=True)
            subprocess.check_call('./initalize_Reference_GW_tilling.sh', shell=True)
            showinfo(
                title='Progressing',
                message='ToOpy is processing the alert (in diagnostic mode) please check command line and output folders')


    def input_gw_event_GW_targeted(self):
            answer = askokcancel(
                title='Confirmation',
                message='Crossmatching alert with galaxy catalog and checking observability from Roque de los Muchachos.',
                icon=WARNING)
            if answer:
                subprocess.check_call('chmod u+r+x initalize_Reference_GW_targeted.sh', shell=True)
                subprocess.check_call('./initalize_Reference_GW_targeted.sh', shell=True)
                showinfo(
                    title='Progressing',
                    message='ToOpy is processing the alert please check command line and output folders')

    def input_gw_event_GW170817_STMOC(self):
            answer = askokcancel(
                title='Confirmation',
                message='Crossmatching alert with galaxy catalog and checking observability from Roque de los Muchachos.',
                icon=WARNING)
            if answer:
                subprocess.check_call('chmod u+r+x initalize_Reference_GW170817_STMOC.sh', shell=True)
                subprocess.check_call('./initalize_Reference_GW170817_STMOC.sh', shell=True)
                showinfo(
                    title='Progressing',
                    message='ToOpy is processing the alert please check command line and output folders')





    def initalize_Reference_GW(self):
        print('x')
        subprocess.check_call('chmod u+r+x initalize_Reference_GW_A.sh', shell=True)
        subprocess.check_call('./initalize_Reference_GW_A.sh', shell=True)
        print('at the end')
        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./documentation_images")
        self.large_test_image = customtkinter.CTkImage(Image.open(os.path.join(image_path, "ToOpy.png")), size=(1000, 300))
        self.home_frame_large_image_label = customtkinter.CTkLabel(self, text="", image=self.large_test_image)
        self.home_frame_large_image_label.grid(row=2, column=0, padx=(20, 0), pady=(20, 0), sticky="nsew")
    '''
    def initalize_Reference_IC170922A_ROI(self):
        subprocess.check_call('chmod u+r+x initalize_Reference_IC170922A_ROI.sh', shell=True)
        subprocess.check_call('./initalize_Reference_IC170922A_ROI.sh', shell=True)
        image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./documentation_images")
        self.large_test_image = customtkinter.CTkImage(Image.open(os.path.join(image_path, "track.png")), size=(1000, 300))
        self.home_frame_large_image_label = customtkinter.CTkLabel(self, text="", image=self.large_test_image)
        self.home_frame_large_image_label.grid(row=2, column=0, padx=(20, 0), pady=(20, 0), sticky="nsew")
        #file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "./")
        #for i in range(0,2):
        #   print(i)
        
    def initalize_Reference_IC170922A_HE(self):
        subprocess.check_call('chmod u+r+x initalize_Reference_IC170922A_HE.sh', shell=True)
        subprocess.check_call('./initalize_Reference_IC170922A_HE.sh', shell=True)

    def initalize_Reference_IC220303A_HE(self):
        subprocess.check_call('chmod u+r+x initalize_Reference_IC220303A_HE.sh', shell=True)
        subprocess.check_call('./initalize_Reference_IC220303A_HE.sh', shell=True)

    def initalize_Reference_IC220303A_ROI(self):
        subprocess.check_call('chmod u+r+x initalize_Reference_IC220303A_ROI.sh', shell=True)
        subprocess.check_call('./initalize_Reference_IC220303A_ROI.sh', shell=True)
    '''



    #def input_gw_event(self):
    #    dialog = customtkinter.CTkInputDialog(text="Add conditions in the following format: ...", title="Confirmation")
    #    print("This is fed into .sh Scripts:", dialog.get_input())

    #def input_neutrino_event_IC170922A_HE(self):
    #    dialog = customtkinter.CTkInputDialog(text="Carefull this will take time to run", title="Progress")
    #    print("This is fed into .sh Scripts:", dialog.get_input())
    #    subprocess.check_call('chmod u+r+x initalize_Reference_IC170922A_HE.sh', shell=True)
    #    subprocess.check_call('./initalize_Reference_IC170922A_HE.sh', shell=True)

    def input_neutrino_event_IC170922A_ROI(self):
        answer = askokcancel(
            title='Confirmation',
            message='Crossmatching alert with galaxy catalog and checking observability from Roque de los Muchachos.',
            icon=WARNING)
        if answer:
            subprocess.check_call('chmod u+r+x initalize_Reference_IC170922A_ROI.sh', shell=True)
            subprocess.check_call('./initalize_Reference_IC170922A_ROI.sh', shell=True)
            showinfo(
                title='Progressing',
                message='ToOpy is processing the alert please check command line and output folders')

    def input_neutrino_event_IC170922A_HE(self):
        answer = askokcancel(
            title='Confirmation',
            message='Running Fermi-LAT Analysis will take some time.',
            icon=WARNING)
        if answer:
            subprocess.check_call('chmod u+r+x initalize_Reference_IC170922A_HE.sh', shell=True)
            subprocess.check_call('./initalize_Reference_IC170922A_HE.sh', shell=True)
            showinfo(
                title='Progressing',
                message='ToOpy is processing the alert please check command line and output folders')

    def input_neutrino_event_IC220303A_ROI(self):
        answer = askokcancel(
            title='Confirmation',
            message='Crossmatching alert with galaxy catalog and checking observability from Roque de los Muchachos.',
            icon=WARNING)
        if answer:
            subprocess.check_call('chmod u+r+x initalize_Reference_IC220303A_ROI.sh', shell=True)
            subprocess.check_call('./initalize_Reference_IC220303A_ROI.sh', shell=True)
            showinfo(
                title='Progressing',
                message='ToOpy is processing the alert please check command line and output folders')

    def input_neutrino_event_IC220303A_HE(self):
        answer = askokcancel(
            title='Confirmation',
            message='Running Fermi-LAT Analysis will take some time (extra time for light curve).',
            icon=WARNING)
        if answer:
            subprocess.check_call('chmod u+r+x initalize_Reference_IC220303A_HE.sh', shell=True)
            subprocess.check_call('./initalize_Reference_IC220303A_HE.sh', shell=True)
            showinfo(
                title='Progressing',
                message='ToOpy is processing the alert please check command line and output folders')
    #def input_neutrino_event_Track(self):
    #    dialog = customtkinter.CTkInputDialog(text="VarInd or HE?", title="Progress")
    #    print("This is fed into .sh Scripts:", dialog.get_input())
    #    if (print(dialog.get_input()) == "VarInd"):
    #    	print('x')

    def input_event_GBM_targeted(self):
        answer = askokcancel(
            title='Confirmation',
            message='Crossmatching alert with galaxy catalog and checking observability from Roque de los Muchachos.',
            icon=WARNING)
        if answer:
            subprocess.check_call('chmod u+r+x initalize_Reference_GBM_targeted.sh', shell=True)
            subprocess.check_call('./initalize_Reference_GBM_targeted.sh', shell=True)
            showinfo(
                title='Progressing',
                message='ToOpy is processing the alert please check command line and output folders')


    def input_event_GBM_STMOC(self):
            answer = askokcancel(
                title='Confirmation',
                message='Initalize STMOC method and check for spatio-temporal coincidences with other alerts.',
                icon=WARNING)
            if answer:
                subprocess.check_call('chmod u+r+x initalize_Reference_GBM_STMOC.sh', shell=True)
                subprocess.check_call('./initalize_Reference_GBM_STMOC.sh', shell=True)
                showinfo(
                    title='Progressing',
                    message='ToOpy is processing the alert please check command line and output folders')



    def input_event_BAT_targeted(self):
        answer = askokcancel(
            title='Confirmation',
            message='Crossmatching alert with galaxy catalog and checking observability from Roque de los Muchachos.',
            icon=WARNING)
        if answer:
            subprocess.check_call('chmod u+r+x initalize_Reference_BAT_targeted.sh', shell=True)
            subprocess.check_call('./initalize_Reference_BAT_targeted.sh', shell=True)
            showinfo(
                title='Progressing',
                message='ToOpy is processing the alert please check command line and output folders')


    def input_event_BAT_STMOC(self):
            answer = askokcancel(
                title='Confirmation',
                message='Initalize STMOC method and check for spatio-temporal coincidences with other alerts.',
                icon=WARNING)
            if answer:
                subprocess.check_call('chmod u+r+x initalize_Reference_BAT_STMOC.sh', shell=True)
                subprocess.check_call('./initalize_Reference_BAT_STMOC.sh', shell=True)
                showinfo(
                    title='Progressing',
                    message='ToOpy is processing the alert please check command line and output folders')

        #argument_list = dialog.get_input()
        #print(argument_list)
        #argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skymap_fits_url), 'VarInd', fermitools_refdata_path, 'no']
        #separator = " "
        #subprocess.check_call("./method_scripts/IceCube_TRACK.sh %s" % separator.join(argument_list), shell=True)

if __name__ == "__main__":
    app = App()
    app.mainloop()





