class TurbidityMonitorRunnable(Runnable):
    def __init__(self,
                 tm: TurbidityMonitor,
                 graph_path: str,
                 camera: Camera,
                 slack_manager,
                 experiment_information: ExperimentInformation,
                 ):
        Runnable.__init__(self, logger=logger)
        self.camera = camera
        self.tm = tm
        self.slack_manager = slack_manager
        self.graph_path = graph_path
        self.experiment_information = experiment_information
        self.time_out_ref: datetime = None
        self.time_out_mins = 10
        self.time_between_measurements = int(60 / tm_n_measurements_per_min)  # seconds between measurements

    def run(self):
        while self.running:
            time.sleep(self.time_between_measurements)  # seconds between measurements

            if pause_monitoring_bool:
                continue

            images = self.camera.take_photos(n=tm_n_images_per_measurement, save_photo=True)
            self.tm.add_measurement(*images)
            self.tm.save_data()
            self.save_current_graph()
            now_time = datetime.now()

            if way_below_dissolved(tm=self.tm):
                slack_manager.post_slack_message('Seems like the turbidity is way below the dissolved reference, '
                                                 'something might be wrong')
                pause_addition_monitoring()
                slack_manager.post_slack_message('If nothing is wrong, you can choose to resume by sending "resume"')

            if pause_addition_bool:
                continue

            # if > time out time since last addition
            if (now_time - self.time_out_ref).seconds > 60 * self.time_out_mins:
                if self.experiment_information.can_add_more_solvent():
                    # maximum solvent hasnt been added in yet, then add more solvent
                    added_solvent = add_solvent(solvent=solvent,
                                                ei=self.experiment_information)
                    if added_solvent is True:
                        try:
                            message = f'{self.time_out_mins} minutes passed, so I will dose in another' \
                                      f' {ei.solvent_addition_volume} mL anyways.\n' \
                                      f'Total volume: {ei.target_volume} mL'
                            slack_manager.post_slack_message(msg=message)
                        except Exception as e:
                            logger.error(e)
                        self.time_out_ref = now_time
                    else:
                        time.sleep(45)
                else:
                    message = f'{self.time_out_mins} minutes passed, but I reached the max volume and cannot add ' \
                              f'any more solvent'
                    slack_manager.post_slack_message(msg=message)
                    self.slack_current_graph()
                    self.slack_current_image()
                    clean_up_vial_fn()
                    time.sleep(1)
                    start_next_vial_fn()
                    time.sleep(60)

            if self.tm.state_changed_to_dissolved() or seems_dissolved(tm=self.tm):
                minutes_to_stir = 5
                minutes_to_check = 5
                message = f'I think the state has changed to dissolved, but to be sure, I will increase the stir rate ' \
                          f'to 900 for {minutes_to_stir} minutes, then I will check the state again, monitoring for {minutes_to_check} minutes'
                slack_manager.post_slack_message(msg=message)
                self.slack_current_graph()
                self.slack_current_image()
                deck.mini_heater_stirrer.target_stir_rate = 900
                time.sleep(minutes_to_stir * 60)
                deck.mini_heater_stirrer.target_stir_rate = self.experiment_information.stir_rate
                # then monitor for to check if it has dissolved
                start_dissolve_check_time = datetime.now()
                current_dissolve_check_time = datetime.now()
                while (current_dissolve_check_time - start_dissolve_check_time).seconds < (minutes_to_check * 60):
                    time.sleep(self.time_between_measurements/2)
                    images = self.camera.take_photos(n=tm_n_images_per_measurement, save_photo=True)
                    self.tm.add_measurement(*images)
                    self.tm.save_data()
                    self.save_current_graph()
                    current_dissolve_check_time = datetime.now()
                if tm.state == TurbidityMonitor.dissolved_state or seems_dissolved(tm=self.tm):
                    message = f'I think the state has really changed to dissolved after 5 minutes of high speed stirring'
                    slack_manager.post_slack_message(msg=message)
                    self.slack_current_graph()
                    self.slack_current_image()
                    self.experiment_information.dissolved = True
                    time.sleep(60)
                else:
                    message = f'I dont think the state is acutally dissolved after all, going to continue monitoring turbidity'
                    slack_manager.post_slack_message(msg=message)
                    self.slack_current_graph()
                    self.slack_current_image()
                    now_time = datetime.now()
            elif tm.state_changed_to_stable():
                message = 'I think the state is steady.\n'
                self.time_out_ref = now_time
                if self.experiment_information.can_add_more_solvent():
                    added_solvent = add_solvent(solvent=solvent,
                                                ei=self.experiment_information)
                    if added_solvent:
                        message += f'I added another {ei.solvent_addition_volume} mL.\n' \
                                   f'Total volume: {ei.target_volume} mL'
                        slack_manager.post_slack_message(msg=message)
                    else:
                        time.sleep(45)
                else:
                    message = 'But I reached the max volume and cannot add any more solvent'
                    slack_manager.post_slack_message(msg=message)
                    self.slack_current_graph()
                    self.slack_current_image()
                    clean_up_vial_fn()
                    time.sleep(1)
                    start_next_vial_fn()
                    time.sleep(60)
            elif tm.state_changed_to_unstable():
                pass

    def start_background_monitoring(self):
        self.time_out_ref = datetime.now()
        self.start()

    def stop_background_monitoring(self):
        self.stop()

    def slack_current_image(self):
        try:
            last_image = camera.last_frame
            last_image_path = str(camera.save_folder.joinpath('last_image.jpg'))
            cv2.imwrite(last_image_path, last_image)
            time.sleep(1)
            self.slack_manager.post_slack_file(filepath=last_image_path,
                                               title='The last image taken',
                                               comment='The last image taken',
                                               )
        except Exception as e:
            print(e)

    def save_current_graph(self):
        try:
            graph_with_liquid_addition(turbidity_data_path=Path(self.tm._turbidity_monitor_data_json_save_path),
                                       liquid_addition_path=solvent_addition_data.csv_path,
                                       save_path=Path(self.graph_path),
                                       tm_parameters=tm_parameters,
                                       )
        except Exception as e:
            print('error saving current graph')
            print(e)

    def slack_current_graph(self):
        try:
            self.save_current_graph()
            time.sleep(1)
            slack_manager.post_slack_file(filepath=self.graph_path,
                                          title='Turbidity vs. time',
                                          comment='Turbidity vs. time graph, the green regions are stable '
                                                  'regions'
                                          )
        except Exception as e:
            print(e)