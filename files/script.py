import tkinter as tk
from tkinter import scrolledtext
import pyautogui
import random
import time
import threading

class HumanTypingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Natural Code Typer")
        self.root.geometry("600x500")
        
        # UI Elements
        self.label = tk.Label(root, text="Paste your code below:", font=("Arial", 11))
        self.label.pack(pady=5)
        
        self.text_area = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=70, height=20)
        self.text_area.pack(pady=10, padx=10, fill=tk.BOTH, expand=True)
        
        self.start_button = tk.Button(root, text="Start Typing (5s Delay)", command=self.start_typing_thread, bg="#4CAF50", fg="white", font=("Arial", 11, "bold"))
        self.start_button.pack(pady=10)
        
        self.status_label = tk.Label(root, text="Status: Ready", fg="gray")
        self.status_label.pack(pady=5)

    def start_typing_thread(self):
        # Run in a separate thread so the GUI doesn't freeze
        typing_thread = threading.Thread(target=self.simulate_typing)
        typing_thread.daemon = True
        typing_thread.start()

    def simulate_typing(self):
        code = self.text_area.get("1.0", tk.END).strip()
        if not code:
            self.status_label.config(text="Status: No code provided!", fg="red")
            return

        # Countdown delay to allow switching to VS Code
        for i in range(5, 0, -1):
            self.status_label.config(text=f"Status: Starting in {i} seconds... Switch to VS Code!", fg="blue")
            self.root.update()
            time.sleep(1)

        self.status_label.config(text="Status: Typing...", fg="green")
        
        i = 0
        while i < len(code):
            char = code[i]
            
            # 1. Simulate a random mistake (approx 2% chance, only on letters/numbers)
            if char.isalnum() and random.random() < 0.02:
                wrong_char = random.choice("abcdefghijklmnopqrstuvwxyz")
                pyautogui.write(wrong_char)
                time.sleep(random.uniform(0.1, 0.25)) # Pause after making mistake
                
                # Backspace to fix it
                pyautogui.press('backspace')
                time.sleep(random.uniform(0.1, 0.2))
            
            # 2. Type the correct character
            pyautogui.write(char)
            
            # 3. Micro-pauses between keystrokes (varies dynamically)
            base_delay = random.uniform(0.05, 0.15)
            
            # Longer pause for punctuation, line breaks, or brackets
            if char in ['.', ',', ';', '\n', '{', '}']:
                base_delay += random.uniform(0.3, 0.6)
                
            time.sleep(base_delay)
            
            # 4. Simulate taking a short "thinking break" (approx 1% chance per character)
            if random.random() < 0.01:
                # Take a break between 1 to 3 seconds
                time.sleep(random.uniform(1.0, 3.0))
                
            i += 1

        self.status_label.config(text="Status: Finished Typing!", fg="gray")

if __name__ == "__main__":
    root = tk.Tk()
    app = HumanTypingApp(root)
    root.mainloop()