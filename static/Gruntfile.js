'use strict';

module.exports = function( grunt ) {

  // Project configuration.
  grunt.initConfig( {

    watch: {
      scss: {
        files: [ 'scss/**/*.scss' ],
        tasks: [ 'compass' ],
        options: {
          livereload: true,
          interrupt: true,
          spawn: true,
        },
      },
      css: {
        files: [ 'css/*.css' ],
        options: {
          livereload: true,
          interrupt: true,
          spawn: true,
        },
      },
    },

    compass: {
      options: {
        sassDir: 'scss',
        cssDir: 'css',
        fontsDir: 'fonts',
        imagesDir: 'images',
        generatedImagesDir: 'images/generated',
        relativeAssets: true,
        noLineComments: true,
        assetCacheBuster: true,
        watch: false,
        require: [ 'breakpoint' ]
      },
      development: {
        options: {
          outputStyle: 'compressed', //nested, expanded, compact, compressed
          environment: 'development',
        }
      },
    },

  });

  grunt.loadNpmTasks( 'grunt-contrib-watch' );
  grunt.loadNpmTasks( 'grunt-contrib-compass' );

  // Tasks.
  grunt.registerTask( 'development', [
    'compass:development',
    'watch',
  ]);

};
