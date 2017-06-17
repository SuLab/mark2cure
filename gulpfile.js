var gulp = require('gulp');
var glibs = require('gulp-load-plugins')({
  pattern: ['gulp-*', 'gulp.*', 'main-bower-files'],
  camelize: true,
});

var path = require('path');
var tinylr = require('tiny-lr');
var server = tinylr();

var livereload = require('gulp-livereload');
var compass = require('gulp-compass');
var csso = require('gulp-csso');
var cleanCSS = require('gulp-clean-css');
var rename = require("gulp-rename");

var SegfaultHandler = require('segfault-handler');
SegfaultHandler.registerHandler();

gulp.task('watch', function () {
  server.listen(35729, function (err) {
    if (err) {
      return console.log(err);
    }
    gulp.watch('static/scss/**/*.scss', ['css']);
    gulp.watch('static/js-src/**/*.js', ['js']);
  });
});

gulp.task('fonts', function() {
  return gulp
    .src(path.join(__dirname, 'node_modules', 'font-awesome', 'fonts/*'))
    .pipe(gulp.dest('static/fonts/'))
    .pipe( livereload(server));
});

gulp.task('js', function() {
  var project_js_files = [
    './node_modules/jquery/dist/jquery.js',
    './static/js-src/libs/jquery-ui.js',
    './node_modules/underscore/underscore.js',
    './node_modules/underscore.string/dist/underscore.string.min.js',

    './node_modules/backbone/backbone.js',
    './node_modules/backbone.radio/build/backbone.radio.js',
    './node_modules/backbone-relational/backbone-relational.js',
    './node_modules/backbone.marionette/lib/backbone.marionette.js',

    './node_modules/tether/dist/js/tether.js',
    './node_modules/bootstrap/dist/js/bootstrap.js',

    './node_modules/d3/d3.js',
    './node_modules/intro.js/intro.js',
    './node_modules/moment/moment.js',
    './node_modules/odometer/odometer.js',
    './node_modules/sifter/sifter.js',
    './node_modules/@rq/sly-scrolling/dist/sly.js',

    './node_modules/sigma/build/sigma.min.js',
    './node_modules/sigma/build/plugins/sigma.parsers.json.min.js',
    './node_modules/sigma/build/plugins/sigma.plugins.filter.min.js',

    './node_modules/raven-js/dist/raven.js',

    './static/js-src/libs/ypet.js',
    './static/js-src/libs/tree.js',
    './static/js-src/tasks/relation-data-var-assigned.js',
    './static/js-src/libs/leaderboard.js',
    './static/js-src/libs/dashboard.js',
    './static/js-src/libs/homepage.js',

    './static/js-src/pages/cloud.js',
    './static/js-src/pages/group_home.js',

    './static/js-src/tasks/relation-synopsis.js',
    './static/js-src/tasks/relation.js',

    './static/js-src/training/entity-recognition/basic.js',
    './static/js-src/training/relation/relation-training.js',
    './static/js-src/training/relation/relation-1.js',

    './static/js-src/demo.js',
    './static/js-src/app.js'
  ];

  return gulp.src(project_js_files)
    // .pipe(glibs.uglify())
    .pipe(glibs.concat('mark2cure.min.js'))
    .pipe(gulp.dest('./static/js/'))
    .on('error', function(err) {
      console.error('Error in compress task', err.toString());
    });
});

gulp.task('css', ['fonts'], function() {
    return gulp.src('./static/scss/style.scss')
      .pipe(compass({
        project: path.join(__dirname, 'static'),
        sass: 'scss',
        import_path: [
          path.join(__dirname, 'node_modules', 'bootstrap', 'scss'),
          path.join(__dirname, 'node_modules', 'font-awesome', 'scss'),
        ],
        time: true,
        debug: true
      }))
      .pipe(csso())
      .pipe(livereload(server))
      // .on('error', sass.logError)
      .pipe(cleanCSS({compatibility: 'ie8'}))
      .pipe(rename({suffix: '.min'}))
      .pipe(gulp.dest('css/'));
});

gulp.task('default', ['js']);

