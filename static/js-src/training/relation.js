//
// Models + Collections
//

RETutorialStep = Backbone.RelationalModel.extend({
  defaults: {
    'order': 0, 'name': '', 'description': '',
    'completed': false
  }
});
RETutorialStepList = Backbone.Collection.extend({
  model: RETutorialStep
});

RETutorial = Backbone.RelationalModel.extend({
  defaults: {
    'order': 0, 'name': '', 'description': '',
    'completed': false,
  },
  relations: [{
    type: 'HasMany',
    key: 'steps',
    relatedModel: RETutorialStep,
    collectionType: RETutorialStepList,
  }]
});
RETutorialList = Backbone.Collection.extend({
  model: RETutorial,
  comparator: 'order'
});

var res = [
  {
    'order': 0,
    'name': 'Using the Interface',
    'description': '',
    'completed': false,
    'steps': [
      {
        'order': 0,
        'name': 'Introduction',
        'description': '',
        'completed': false,
        'paragraphs': [
          {
            'order': 0,
            'text': 'There are many different kinds of concepts in biomedical text. With the concept recognition task, you identify and tag concepts like genes, disease, and treatments in biomedical abstracts. In this module, you learn how to use the relationship extraction tool to annotate the relationship between different concepts in biomedical text and earn your Relationship Marking Skill.',
            'muted': false
          }
        ]
      }
    ]
  },
  {
    'order': 1,
    'name': 'Rules for Relationship Extraction',
    'description': '',
    'completed': false,
    'steps': [
    ]
  },

  {
    'order': 2,
    'name': 'Learning Relationships',
    'description': '',
    'completed': false,
    'steps': [
      {
        'order': 1,
        'description': 'Genes',
        'completed': false,
        'steps': []
      },

      {
        'order': 0,
        'description': 'Diseases',
        'completed': false,
        'steps': []
      },

      {
        'order': 2,
        'description': 'Treatments',
        'completed': false,
        'steps': []
      },

      {
        'order': 3,
        'description': 'Multi-Marking',
        'completed': false,
        'steps': []
      }
    ]
  },
];

//
// Views
//

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

  regions: {
    'interactive': '#interactive'
  },

  events: {
    // 'click @ui.next': function() { this.triggerMethod('step:complete:2', this); },
  },

  onRender: function() {
    // var main = this.getRegion('interactive');
    //
    // var collection = new REExtractionList([{
    //     "id": 0,
    //     "document": 0,
    //     "relation_type": "g_d",
    //     "concepts": {
    //       "c2": {
    //         "text": "Physics",
    //         "type": "d",
    //         "id": "0"
    //       },
    //       "c1": {
    //         "text": "Volunteers",
    //         "type": "g",
    //         "id": "0"
    //       }
    //     },
    //     "user_completed": false}]);
    // var tree = new Tree({'collection': collection});
    // main.show(tree.render());

    // var sub_view = new Tree({'collection': collection});
    // sub_view.render();

  //
  //   Tree['start'].onShow = function(m, c) {
  //     c.$el.css({'opacity': .25});
  //     c.$el.append('<div id="tree-cover"></div>');
  //
  //     #<{(| (TODO) Integrate into Tree app |)}>#
  //     var current_relationship = collection.findWhere({'current': true});
  //     add_relation_classes(current_relationship);
  //
  //     $el = $('#c1');
  //     $el.popover('hide');
  //     $el.popover({
  //       container: 'body',
  //       html: true,
  //       animation: false,
  //       content: function() {
  //           return 'A concept is a category to which a term can be assigned';
  //       },
  //       placement: 'top'
  //     });
  //     $el.popover('show');
  //
  //     $el2 = $('#c2');
  //     $el2.popover('hide');
  //     $el2.popover({
  //       container: 'body',
  //       html: true,
  //       animation: false,
  //       content: function() {
  //           return 'Eg- ‘Physics’ is part of the ‘field of science’ concept';
  //       },
  //       placement: 'top'
  //     });
  //     $el2.popover('show');
  //   }

  }
});

RETrainingNavigationProgressItem = Backbone.Marionette.View.extend({
  template: _.template('<%= name %>'),
  tagName: 'li',
  className: 'breadcrumb-item'
});

RETrainingNavigationProgress = Backbone.Marionette.CollectionView.extend({
  childView: RETrainingNavigationProgressItem,
  tagName: 'ol',
  className: 'breadcrumb'
});

RETrainingNavigation = Backbone.Marionette.View.extend({
  template: '#re-training-navigation',

  regions: {
    'steps': '#progress-steps'
  },

  onRender: function() {
    this.showChildView('steps', new RETrainingNavigationProgress({collection: this.collection}));
  }
});


RETrainingText = Backbone.Marionette.View.extend({
  template: '#re-training-text',
});


RETrainingAction = Backbone.Marionette.View.extend({
  template: '#re-training-action',
});


RETrainingFooter = Backbone.Marionette.View.extend({
  template: '#re-training-footer',
});


RETraining = Backbone.Marionette.View.extend({
  template: '#re-training-template',

  regions: {
    'navigation': '#re-training-navigation',
    'text': '#re-training-text',
    'action': '#re-training-action',
    'footer': '#re-training-footer',
  },

  // childViewEvents: {
  //   'step:complete:1': function(childView) {
  //     this.emptyRegions();
  //     this.showChildView('two', new Step1());
  //   },
  //   'step:complete:2': function(childView) {
  //     this.emptyRegions();
  //     this.showChildView('three', new Step2());
  //   },
  // },
  //

  initialize: function() {
    this.collection = new RETutorialList(res);
    console.log(this.collection);
  },

  onRender: function() {
    this.showChildView('navigation', new RETrainingNavigation({collection: this.collection}));
    this.showChildView('text', new RETrainingText({collection: this.collection}));
    this.showChildView('action', new RETrainingAction({collection: this.collection}));
    this.showChildView('footer', new RETrainingFooter({collection: this.collection}));
  }

});

