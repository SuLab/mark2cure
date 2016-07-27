var gulp = require('gulp');
var compass = require('gulp-compass');
var gulpif = require('gulp-if');
var livereload = require('gulp-livereload');
var cleanCSS = require('gulp-clean-css');
var gutil = require('gulp-util');
var csso = require('gulp-csso');
var sass = require('gulp-sass');
var rename = require('gulp-rename');

var path = require('path');
var tinylr = require('tiny-lr');
var server = tinylr();

var SegfaultHandler = require('segfault-handler');
SegfaultHandler.registerHandler();

var PRODUCTION_MODE = gutil.env.production;


gulp.task('watch', function () {
  server.listen(35729, function (err) {
    if (err) {
      return console.log(err);
    }
    gulp.watch('static/scss/**/*.scss', ['css']);
  });
});

gulp.task('fonts', function() {
  return gulp
    .src(path.join(__dirname, 'bower_components', 'font-awesome', 'fonts/*'))
    .pipe(gulp.dest('static/fonts/'))
    .pipe( livereload(server));
});

gulp.task('css', ['fonts'], function() {
    return gulp.src('./static/scss/style.scss')
      .pipe(compass({
        project: path.join(__dirname, 'static'),
        sass: 'scss',
        import_path: [
          path.join(__dirname, 'bower_components', 'bootstrap', 'scss'),
          path.join(__dirname, 'bower_components', 'font-awesome', 'scss'),
        ],
        time: true,
        debug: true
      }))
      .pipe(csso())
      .pipe(livereload(server))
      .on('error', sass.logError)
      .pipe(cleanCSS({compatibility: 'ie8'}))
      .pipe(rename({suffix: '.min'}))
      .pipe(gulp.dest('css/'));
});

gulp.task('default', ['css']);

