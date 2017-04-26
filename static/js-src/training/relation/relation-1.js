IntroStep = Backbone.Marionette.View.extend({
  // The initial page, classify + describe the activity
  template: '#intro-step',

  ui: {
    'next': '#next',
  },

  events: {
    'click @ui.next': function() { this.triggerMethod('step:complete:1', this); },
  },

});

Step1 = Backbone.Marionette.View.extend({
  // The initial page, classify + describe the activity
  template: '#step-one',

  ui: {
    'next': '#next',
  },

  events: {
    'click @ui.next': function() { this.triggerMethod('step:complete:2', this); },
  },

  onRender: function() {
    relation_task_settings = {};
    relation_task_settings['data'] = relation_data;
    relation_task_settings['fade'] = true;

    task_data = [{
      "id": 0,
      "document": 0,
      "relation_type": "g_d",
      "concepts": {
        "c2": {
          "text": "Physics",
          "type": "d",
          "id": "0"
        },
        "c1": {
          "text": "Volunteers",
          "type": "g",
          "id": "0"
        }
      },
      "user_completed": false
    }];


    var RelationApp = Backbone.Marionette.Application.extend({
      regions: {
        start: '#tree-insert'
      },

      onStart: function() {
        Backbone.Radio.DEBUG = true;
        RelationApp['convoChannel'] = Backbone.Radio.channel('convo');

        var main = this.getRegion('start');  // Has all the properties of a `Region`
        main.show(new SomeView());
      }
    });
    var app = new RelationApp();

    // Tree['start'].onShow = function(m, c) {
    //   if(relation_task_settings['fade']) {
    //     c.$el.css({'opacity': .25});
    //     c.$el.append('<div id="tree-cover"></div>');
    //   }

    /* Initalize the page by using the task_data */
    var collection;
    /* Initalize the Application, but load events later */
    app.start();

    /* Onload request all relation tasks to complete */
    collection = new RelationTaskCollection(task_data);
    /* Sort the collection by C1 on initial data load */
    var new_collection = _.sortBy(collection.models, function(c) {
      return c.attributes.concepts.c1.text;
    });
    collection.models = new_collection;
    collection.next();

    var current_relationship = collection.findWhere({'current': true});
    add_relation_classes(current_relationship);


    $el = $('#c1');
    $el.popover('hide');
    $el.popover({
      container: 'body',
      html: true,
      animation: false,
      content: function() {
          return 'A concept is a category to which a term can be assigned';
      },
      placement: 'top'
    });
    $el.popover('show');

    $el2 = $('#c2');
    $el2.popover('hide');
    $el2.popover({
      container: 'body',
      html: true,
      animation: false,
      content: function() {
          return 'Eg- ‘Physics’ is part of the ‘field of science’ concept';
      },
      placement: 'top'
    });
    $el2.popover('show');

  }
});



RelationTrainingOne = Backbone.Marionette.View.extend({
  el: '#relation-training-one',
  template: false,

  regions: {
    'one': '#step-1',
    'two': '#step-2',
    'three': '#step-3',
    'four': '#step-4',
    'five': '#step-5'
  },

  childViewEvents: {
    'step:complete:1': function(childView) {
      this.emptyRegions();
      this.showChildView('two', new Step1());
    },
    'step:complete:2': function(childView) {
      this.emptyRegions();
      this.showChildView('three', new Step2());
    },
  },

  onRender: function() {
    // Initial load of the page
    this.showChildView('one', new IntroStep());
  }

});

if($('#relation-training-one').length) {
  var sub_view = new RelationTrainingOne();
  sub_view.render();
}
