#include <iostream>
#include <dirent.h>
#include <sys/stat.h>

/**
 * @brief Create a directory if it doesn not exist
 * 
 * @param dirname 
 */
void 
create_directory(const char *dirname)
{
    if(opendir(dirname) == NULL){
        mkdir(dirname,0777);
    } 
}

/**
 * @brief print progress bar to the screen
 * @param percentage percentage of current progress
*/
void print_progressbar(float percentage) {
    // Ensure percentage is within [0.0, 100.0]
    percentage = (percentage < 0.0) ? 0.0 : (percentage > 100.0) ? 100.0 : percentage;

    // Determine the number of characters to represent the progress bar
    int numChars = (int)(percentage / 2.0);

    // Print the progress bar
    printf("[");
    for (int i = 0; i < 50; ++i) {
        if (i < numChars) {
            printf("=");
        } else if (i == numChars) {
            printf(">");
        } else {
            printf(" ");
        }
    }
    printf("] %.1f%%\r", percentage);  // Use carriage return to overwrite the line

    // Flush the output to ensure immediate display
    fflush(stdout);
}