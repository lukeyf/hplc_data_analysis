import os
import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


class Watcher:
    DIRECTORY_TO_WATCH = r"D:\Chemstation\1\Data"
    def __init__(self):
        self.observer = Observer()

    def run(self):
        event_handler = Handler()
        self.observer.schedule(event_handler, self.DIRECTORY_TO_WATCH, recursive=True)
        self.observer.start()
        i = 3000
        try:
            while deck.reaction_folder=='' and i >= 1:
                time.sleep(1)
                i -=1
        except:
            self.observer.stop()
            print("Error")
        self.observer.stop()



class Handler(FileSystemEventHandler):
    @staticmethod
    def on_any_event(event):
        if event.is_directory:
            return None

        elif event.event_type == 'created':
            try:

                deck.reaction_folder = take_n_level_dir(event.src_path,5)
                print(deck.reaction_folder)
            except:
                print('pass')

        elif event.event_type == 'modified':
            # Taken any action here when a file is modified.
            pass

def take_n_level_dir(path:str,n:int):
    """
    works for mac system
    :param path:
    :param n:
    :return:
    """
    split = path.split('\\')
    x = split[n+1]
    join = '\\'.join(split[:n+1])
    return join


if __name__ == '__main__':
    w = Watcher()
    w.run()