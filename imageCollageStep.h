#pragma once

#include <sys/stat.h>
#include <vector>
#include <algorithm>
#include <cstdlib> // для std::system
#include <filesystem>

#include "source/manageClasses/baseStep.h"
#include "MM_project_utils.h"
#include "MM_project_const.h"

class ImageCollageStep : public BaseStep {
public:
    /**
     * @param stepName Имя шага
     * @param inputDir Директория с исходными PNG
     * @param outputDir Директория для сохранения коллажей
     */
    ImageCollageStep(const std::string& stepName,
                     const std::string& inputDir,
                     const std::string& outputDir,
                     const std::string& archiveName = "")
        : BaseStep(stepName, {}, ""), 
          inputDir(inputDir), 
          outputDir(outputDir),
          archiveName(archiveName)
    {
        /*
        // Для работы автоматики isDone() назначим "контрольный файл" - 
        // последний файл, который должен появиться в результате работы.
        // Если он есть, BaseStep решит, что шаг выполнен.
        std::string lastCollage = outputDir + "summary_q2_" + std::to_string(NUMBER_Q2) + 
                                  "_w_" + std::to_string(NUMBER_W) + 
                                  "_ct_" + std::to_string(NUMBER_COS_THETA) + ".png";
        setOutputFile(lastCollage);
        */
    }

protected:
    std::string inputDir;
    std::string outputDir;
    std::string archiveName;

    bool initialize() override {
        /*
        // Проверяем, заканчивается ли путь на слеш, если нет - добавляем
        if (!inputDir.empty() && inputDir.back() != '/') inputDir += '/';
        if (!outputDir.empty() && outputDir.back() != '/') outputDir += '/';

        // Создаем папку для коллажей
        // 0777 — права доступа (rwxrwxrwx)
        mkdir(outputDir.c_str(), 0777);
        */
        
        return true;
    }

    bool run() override {
        log("Starting collage generation...", LogLevel::Info);

        for (int iq = 0; iq < NUMBER_Q2; ++iq) {
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    
                    std::string outputName = outputDir + "summary_q2_" + std::to_string(iq+1) + 
                                            "_w_" + std::to_string(iw+1) + 
                                            "_ct_" + std::to_string(ict+1) + ".png";

                    std::string listFileName = outputDir + "temp_file_list.txt";
                    std::ofstream listFile(listFileName);
                    
                    if (!listFile.is_open()) {
                        log("Could not create temporary list file!", LogLevel::Error);
                        return false;
                    }

                    int foundFilesCount = 0; 

                    for (int iphi = 0; iphi < NUMBER_PHI; ++iphi) {
                        std::string fileName = inputDir + makeCellName(iq, iw, ict, iphi) + ".png";
                        
                        if (std::filesystem::exists(fileName)) {
                            // Файл есть — пишем его путь (в кавычках)
                            listFile << "\"" << fileName << "\"\n";
                            foundFilesCount++;
                        } else {
                            // ФАЙЛА НЕТ — вставляем пустышку для сохранения позиции
                            // null: — специальное ключевое слово ImageMagick
                            listFile << "null:\n"; 
                        }
                    }
                    listFile.close();

                    // Если нашли хотя бы одно изображение — запускаем монтаж
                    if (foundFilesCount > 0) {
                        log("Creating collage " + outputName + " (" + std::to_string(foundFilesCount) + " images)", LogLevel::Debug);
                        
                        std::string command = "montage @" + listFileName + " -tile 5x2 -geometry +2+2 " + outputName;
                        int result = std::system(command.c_str());
                        
                        if (result != 0) {
                            log("Montage failed for: " + outputName, LogLevel::Warning);
                        }
                    } else {
                        log("No images found for bin Q2:" + std::to_string(iq+1) + 
                            " W:" + std::to_string(iw+1) + " CT:" + std::to_string(ict+1) + ". Skipping.", LogLevel::Info);
                    }

                    std::remove(listFileName.c_str());
                }
            }
        }
        return true;
    }

    bool finalize() override {
        log("All collages processed.", LogLevel::Info);

        if (!archiveName.empty()) {
            log("Creating archive: " + archiveName, LogLevel::Info);

            // Оборачиваем пути в кавычки на случай пробелов в именах папок или архива
            // Было: ... + archiveName + " -C " + outputDir + " ."
            // Стало: ... + "\"" + archiveName + "\" -C \"" + outputDir + "\" ."
            
            std::string tarCmd = "tar -czf \"" + archiveName + 
                                "\" -C \"" + outputDir + "\" .";
            
            log("Executing: " + tarCmd, LogLevel::Debug);
            
            int result = std::system(tarCmd.c_str());
            if (result == 0) {
                log("Archive created successfully: " + archiveName, LogLevel::Info);
            } else {
                log("Failed to create archive. Check if 'tar' is installed and paths are correct.", LogLevel::Error);
                return false;
            }
        }

        return true;
    }

    /*
    // Переопределяем isValid, так как входные файлы у нас — это динамическая папка
    bool isValid() const override {
        struct stat info;
        if (stat(inputDir.c_str(), &info) != 0) {
            log("Input directory does not exist: " + inputDir, LogLevel::Error);
            return false;
        }
        return (info.st_mode & S_IFDIR);
    }
    */
};